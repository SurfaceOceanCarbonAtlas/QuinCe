package uk.ac.exeter.QuinCe.data.Dataset;

import java.time.LocalDateTime;
import java.time.temporal.ChronoUnit;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.sql.DataSource;

import org.apache.commons.lang3.mutable.MutableInt;

import uk.ac.exeter.QuinCe.data.Files.DataFile;
import uk.ac.exeter.QuinCe.data.Files.DataFileDB;
import uk.ac.exeter.QuinCe.data.Files.DataFileException;
import uk.ac.exeter.QuinCe.data.Files.DataFileLine;
import uk.ac.exeter.QuinCe.data.Instrument.FileDefinition;
import uk.ac.exeter.QuinCe.data.Instrument.FileDefinitionException;
import uk.ac.exeter.QuinCe.data.Instrument.Instrument;
import uk.ac.exeter.QuinCe.data.Instrument.Calibration.CalibrationNotValidException;
import uk.ac.exeter.QuinCe.data.Instrument.Calibration.CalibrationSet;
import uk.ac.exeter.QuinCe.data.Instrument.Calibration.SensorCalibrationDB;
import uk.ac.exeter.QuinCe.data.Instrument.DataFormats.DateTimeSpecificationException;
import uk.ac.exeter.QuinCe.data.Instrument.DataFormats.PositionException;
import uk.ac.exeter.QuinCe.data.Instrument.RunTypes.RunTypeCategory;
import uk.ac.exeter.QuinCe.data.Instrument.SensorDefinition.SensorAssignment;
import uk.ac.exeter.QuinCe.data.Instrument.SensorDefinition.SensorType;
import uk.ac.exeter.QuinCe.utils.DatabaseException;
import uk.ac.exeter.QuinCe.utils.DateTimeUtils;
import uk.ac.exeter.QuinCe.utils.ExtendedMutableInt;
import uk.ac.exeter.QuinCe.utils.MissingParamException;
import uk.ac.exeter.QuinCe.utils.RecordNotFoundException;

/**
 * Class to load the raw data for a data set
 * @author Steve Jones
 *
 */
public abstract class DataSetRawData {

  /**
   * Line index indicating the we've navigated past the end of the file
   */
  protected static final int EOF_VALUE = -1;

  /**
   * Line index indicating that file processing has not started
   */
  protected static final int NOT_STARTED_VALUE = -2;

  /**
   * Line index indicating that we've navigated past the end of the file
   */
  protected static final ExtendedMutableInt EOF = new ExtendedMutableInt(EOF_VALUE);

  /**
   * Line index indicating that we haven't started processing the file
   */
  protected static final ExtendedMutableInt NOT_STARTED = new ExtendedMutableInt(NOT_STARTED_VALUE);

  /**
   * Averaging mode for no averaging
   */
  public static final int AVG_MODE_NONE = 0;

  /**
   * Human readable string for the no-averaging mode
   */
  public static final String AVG_MODE_NONE_NAME = "None";

  /**
   * Averaging mode for averaging every minute
   */
  public static final int AVG_MODE_MINUTE = 1;

  /**
   * Human-readable string for the every-minute averaging mode
   */
  public static final String AVG_MODE_MINUTE_NAME = "Every minute";

  /**
   * The data set to which this data belongs
   */
  private DataSet dataSet;

  /**
   * The instrument to which the data set belongs
   */
  private Instrument instrument;

  /**
   * The file definitions for the data set
   */
  private List<FileDefinition> fileDefinitions;

  /**
   * The data files in the data set
   */
  private List<DataFile> dataFiles;

  /**
   * The data extracted from the data files that are encompassed by the data set
   */
  protected List<List<DataFileLine>> data;

  /**
   * The rows from each file that are currently being processed.
   * These rows are discovered using #nextRecord
   */
  protected List<List<Integer>> selectedRows = null;

  /**
   * Current position pointers for each file definition
   */
  protected List<Integer> linePositions = null;

  /**
   * A calibration set for this DataSet, initialized with an empty
   * calibrationSet
   */
  private CalibrationSet calibrationSet;

  /**
   * Constructor - loads and extracts the data for the data set
   * @param dataSource A data source
   * @param dataSet The data set
   * @param instrument The instrument to which the data set belongs
   * @throws RecordNotFoundException If no data files are found within the data set
   * @throws DatabaseException If a database error occurs
   * @throws MissingParamException If any required parameters are missing
   * @throws DataFileException If the data cannot be extracted from the files
   */
  protected DataSetRawData(DataSource dataSource, DataSet dataSet, Instrument instrument) throws MissingParamException, DatabaseException, RecordNotFoundException, DataFileException {

    this.dataSet = dataSet;
    this.instrument = instrument;

    fileDefinitions = new ArrayList<FileDefinition>();
    data = new ArrayList<List<DataFileLine>>();
    selectedRows = new ArrayList<List<Integer>>();
    linePositions = new ArrayList<Integer>();

    for (FileDefinition fileDefinition : instrument.getFileDefinitions()) {
      fileDefinitions.add(fileDefinition);

      dataFiles = DataFileDB.getFiles(dataSource, fileDefinition, dataSet.getStart(), dataSet.getEnd());
      data.add(extractData(dataFiles));

      selectedRows.add(null);
      linePositions.add(NOT_STARTED_VALUE);
    }
    CalibrationSet calibrations = new SensorCalibrationDB()
        .getMostRecentCalibrations(dataSource, instrument.getDatabaseId(),
            dataSet.getStart());
    if (calibrations.isValid()) {
      setCalibrationSet(calibrations);
    }
    else {
      throw new CalibrationNotValidException("Missing valid calibration");
    }
  }

  /**
   * Extract the rows from the data files that are encompassed by a data set
   * @param dataFiles The data files for the definition
   * @throws DataFileException If the data cannot be extracted from the file
   */
  private List<DataFileLine> extractData(List<DataFile> dataFiles) throws DataFileException {

    // The data from the files
    ArrayList<DataFileLine> fileData = new ArrayList<DataFileLine>();

    // Loop through each file
    for (DataFile file : dataFiles) {

      int currentLine = DataFileException.NO_LINE_NUMBER;
      try {
        // Size the data array
        fileData.ensureCapacity(fileData.size() + file.getRecordCount());

        // Skip any lines before the data set start date
        currentLine = file.getFirstDataLine();
        while (file.getDate(currentLine).isBefore(dataSet.getStart())) {
          currentLine++;
        }

        // Copy lines until either the end of the file or we pass the data set's end date
        while (currentLine < file.getContentLineCount()) {
          if (file.getDate(currentLine).isAfter(dataSet.getEnd())) {
            break;
          } else {
            fileData.add(new DataFileLine(file, currentLine));
            currentLine++;
          }
        }
      } catch (DateTimeSpecificationException e) {
        throw new DataFileException(file.getDatabaseId(), currentLine, e);
      }
    }

    return fileData;
  }

  /**
   * Get the available averaging modes as a map
   * @return The averaging modes
   */
  public static Map<String, Integer> averagingModes() {
    // TODO Move this to a static block
    LinkedHashMap<String, Integer> map = new LinkedHashMap<String, Integer>();

    map.put(AVG_MODE_NONE_NAME, AVG_MODE_NONE);
    map.put(AVG_MODE_MINUTE_NAME, AVG_MODE_MINUTE);

    return map;
  }

  /**
   * Navigate to the next record in the data set.
   *
   * <p>
   *   This will collect all rows that are relevant under the averaging scheme
   *   for the data set.
   * </p>
   *
   * @return {@code true} if a new record is discovered; {@code false} if there are no more records in the data set
   * @throws DataSetException If an error occurs during record selection
   * @throws DataFileException If any data cannot be extracted
   */
  public DataSetRawDataRecord getNextRecord() throws DataSetException, DataFileException {

    DataSetRawDataRecord record = null;

    // If a record is found but errors, then we search for another one.
    // If a record is not found at all, it doesn't matter.
    boolean finishedSearch = false;

    while (!finishedSearch) {
      finishedSearch = true;
      boolean found = false;

      // If there is a previous selection, find the file with
      // the smallest increment to the next row, and start searching from there. This
      // ensures we maximise the number of records from the file.
      //
      // Only records from files with Run Types (i.e. those containing the base fundamental measurements)
      // are checked, otherwise the same fundamental measurement will be used multiple times.
      boolean begin = true;
      int currentFile = 0;

      try {
        long smallestIncrement = Long.MAX_VALUE;
        int smallestIncrementFile = 0;
        for (int i = 0; i < fileDefinitions.size(); i++) {
          if (fileDefinitions.get(i).hasRunTypes()) {
            int currentRow = linePositions.get(i);

            if (currentRow != NOT_STARTED_VALUE && currentRow + 1 < getFileSize(i)) {
              DataFileLine currentLine = getLine(i, currentRow);
              DataFileLine nextLine = getLine(i, currentRow + 1);

              long increment = ChronoUnit.SECONDS.between(currentLine.getDate(), nextLine.getDate());
              if (increment < smallestIncrement) {
                smallestIncrement = increment;
                smallestIncrementFile = i;
              }
            }
          }
        }

        currentFile = smallestIncrementFile;
      } catch (Exception e) {
        throw new DataSetException(e);
      }

      while (begin || !allRowsMatch()) {

        if (!selectNextRow(currentFile)) {
          // We've gone off the end of this file, so we can't select a record
          // across all files
          break;
        } else {
          // If the rows aren't equal, then reset and continue from here
          if (!selectedRowsMatch(currentFile)) {
            resetOtherFiles(currentFile);
          }

          // Select the next file
          currentFile++;
          if (currentFile == fileDefinitions.size()) {
            currentFile = 0;
          }
        }

        begin = false;
      }

      if (allRowsMatch()) {
        found = true;
      }

      if (found) {
        try {
          record = createRecord();
        } catch (DataFileException e) {

          // Creating that record failed. Log the message and try another one.
          System.out.println(e.getMessage());
          finishedSearch = false;
        }
      }
    }

    return record;
  }

  /**
   * Determine whether or not matching lines have been found
   * for all files in the data set
   * @return {@code true} if all files have matching lines; {@code false} otherwise
   */
  private boolean allRowsMatch() throws DataSetException {
    boolean match = true;

    for (int i = 0; match && i < selectedRows.size(); i++) {
      if (null == selectedRows.get(i)) {
        match = false;
      }
    }

    if (match) {
      match = selectedRowsMatch(0);
    }

    return match;
  }

  /**
   * Clear the selected rows from all files except the specified one
   * @param file The source file that should not be cleared
   */
  private void resetOtherFiles(int file) {
    for (int i = 0; i < selectedRows.size(); i++) {
      if (i != file) {
        selectedRows.set(i, null);
      }
    }
  }

  /**
   * Determine whether all the files that have selected rows
   * match the selected rows in the specified file. Files that
   * do not have rows assigned are not checked.
   *
   * @param file The file to compare
   * @return {@code true} if the all the assigned files match; {@code false} otherwise
   * @throws DataSetException If an error occurs during the checks
   */
  private boolean selectedRowsMatch(int file) throws DataSetException {
    boolean match = true;

    for (int i = 0; match && i < selectedRows.size(); i++) {
      if (i != file) {
        match = lineSelectionsMatch(file, i);
      }
    }

    return match;
  }

  /**
   * Compare the selected rows for two files to see if they match
   * @param file1 The first file
   * @param file2 The second file
   * @return {@code true} if the selected rows match; {@code false} if they do not
   * @throws DataSetException If an error occurs during the examination
   */
  protected abstract boolean lineSelectionsMatch(int file1, int file2) throws DataSetException;

  /**
   * Select the next row(s) in the specified file
   * @param fileIndex The file index
   * @return {@code true} if the next row is found; {@code false} if the end of the file has been reached
   * @throws DataSetException If an error occurs during row selection
   */
  protected abstract boolean selectNextRow(int fileIndex) throws DataSetException;

  /**
   * Get a selected line from a file other than the specified file.
   *
   * The file from which the line is taken is not defined. If no other
   * files have selected lines, the method returns {@code null}.
   * @return A selected line from another file
   */
  protected DataFileLine getOtherSelectedLine(int fileIndex) {

    DataFileLine line = null;

    for (int i = 0; i < selectedRows.size(); i++) {
      if (i != fileIndex) {
        List<Integer> lines = selectedRows.get(i);
        if (null != lines) {
          line = getLine(i, lines.get(0));
          break;
        }
      }
    }

    return line;
  }

  /**
   * Reset the line processing back to the start of the files
   */
  public void reset() {
    selectedRows = new ArrayList<List<Integer>>(fileDefinitions.size());
    linePositions = new ArrayList<Integer>(fileDefinitions.size());

    for (int i = 0; i < fileDefinitions.size(); i++) {
      selectedRows.add(null);
      linePositions.add(NOT_STARTED_VALUE);
    }
  }

  /**
   * Find the next usable line in the specified file, starting at the stored line position
   * for that file. See {@link #getNextLine(int, int)} for full details.

   * @param fileIndex The file whose next line is to be located
   * @return The next line index
   * @throws FileDefinitionException If the Run Types for the file are invalid
   * @throws DataFileException If data cannot be extracted from the data lines
   * @see #getNextLine(int, int)
   */
  protected int getNextLine(int fileIndex) throws DataFileException, FileDefinitionException {
    return getNextLine(fileIndex, linePositions.get(fileIndex));
  }

  /**
   * Find the next usable line in the specified file, starting at the specified line.
   *
   * <p>
   *   For files with Run Types, the method skips ignored run types and lines
   *   within the pre- and post- flushing times of the instrument. For files
   *   without Run Types, every line is returned.
   * </p>
   *
   * @param fileIndex The file whose next line is to be located
   * @param startLine The starting point for the search
   * @return The next line index
   * @throws FileDefinitionException If the Run Types for the file are invalid
   * @throws DataFileException If data cannot be extracted from the data lines
   */
  protected int getNextLine(int fileIndex, int startLine) throws DataFileException, FileDefinitionException {

    // TODO This could probably be better written as a state machine

    // TODO Validate date, longitude, latitude for each line. If invalid, record in the data set and move on.

    ExtendedMutableInt result = EOF;
    boolean finished = false;

    int fileLastLine = getFileSize(fileIndex) - 1;

    FileDefinition fileDefinition = fileDefinitions.get(fileIndex);
    ExtendedMutableInt currentLineIndex = new ExtendedMutableInt(startLine);

    // If the file does not have Run Types, return the next line
    if (!fileDefinition.hasRunTypes()) {

      if (currentLineIndex.equals(NOT_STARTED)) {
        result = new ExtendedMutableInt(0);
      } else {
        result = currentLineIndex.incrementedClone();
        if (result.greaterThanOrEqualTo(fileLastLine)) {
          result = EOF;
        }
      }
      finished = true;
    } else {
      // Get the Run Type of the current line
      RunTypeCategory currentRunType = null;
      if (!currentLineIndex.equals(NOT_STARTED)) {
        currentRunType = getLine(fileIndex, currentLineIndex).getRunTypeCategory();
      }

      // Get the next line index
      ExtendedMutableInt nextLineIndex;
      if (currentLineIndex.equals(NOT_STARTED)) {
        nextLineIndex = new ExtendedMutableInt(0);
      } else {
        nextLineIndex = currentLineIndex.incrementedClone();
      }

      if (nextLineIndex.greaterThanOrEqualTo(fileLastLine)) {
        result = EOF;
        finished = true;
      } else {
        String nextRunType = getLine(fileIndex, nextLineIndex).getRunType();

        if (nextRunType.equals(currentRunType)) {
          // Check to see if we're in the post-flushing period
          // If we aren't, we have found the next line
          // If we are, the method will skip us ahead to the next
          // line with a different Run Type, and we can continue from there
          // in the next if block below
          if (!inPostFlushingPeriod(fileIndex, nextRunType, nextLineIndex)) {
            result = nextLineIndex;
            finished = true;
          }
        }

        if (!finished) {

          // Get the new Run Type
          nextRunType = getLine(fileIndex, nextLineIndex).getRunType();

          // We're in a new Run Type. Skip any Ignored lines
          while (!finished && getLine(fileIndex, nextLineIndex).isIgnored()) {
            nextLineIndex.increment();
            if (nextLineIndex.greaterThanOrEqualTo(fileLastLine)) {
              result = EOF;
              finished = true;
            }
          }

          // Assuming we didn't fall off the end of the file...
          if (!finished) {
            String newRunType = getLine(fileIndex, nextLineIndex).getRunType();

            // We are at the start of a new run type. Skip over the pre-flushing stage
            boolean preFlushingTimeProcessed = false;

            while (!preFlushingTimeProcessed) {
              if (skipPreFlushingTime(fileIndex, newRunType, nextLineIndex)) {
                result = nextLineIndex;
                finished = true;
                preFlushingTimeProcessed = true;
              } else {
                while (getLine(fileIndex, nextLineIndex).isIgnored()) {
                  nextLineIndex.increment();
                  if (nextLineIndex.greaterThanOrEqualTo(getFileSize(fileIndex))) {
                    result = EOF;
                    finished = true;
                    preFlushingTimeProcessed = true;
                  }
                }

                newRunType = getLine(fileIndex, nextLineIndex).getRunType();
              }
            }
          }
        }
      }
    }

    return result.getValue();
  }

  /**
   * Determine whether or not the specified line is within the post-flushing period for the instrument
   * @param fileIndex The file
   * @param runType The current Run Type
   * @param line The line being examined
   * @return {@code true} if the line is within the post-flushing period; {@code false} if it is not
   * @throws DataFileException If the line data cannot be processed
   * @throws FileDefinitionException If an invalid Run Type is detected
   */
  private boolean inPostFlushingPeriod(int fileIndex, String runType, ExtendedMutableInt line) throws DataFileException, FileDefinitionException {

    boolean inPostFlushing = false;

    if (instrument.getPostFlushingTime() == 0) {
      inPostFlushing = false;
    } else {

      try {
        LocalDateTime lineDate = getLine(fileIndex, line).getDate();

        boolean finished = false;
        int previousLineIndex = line.getValue();
        LocalDateTime previousDate = lineDate;

        while (!finished) {

          int nextLineIndex = previousLineIndex + 1;
          if (nextLineIndex >= getFileSize(fileIndex)) {
            inPostFlushing = (DateTimeUtils.secondsBetween(lineDate, previousDate) <= instrument.getPostFlushingTime());
            finished = true;
          } else {
            DataFileLine nextLine = getLine(fileIndex, nextLineIndex);
            LocalDateTime nextDate = nextLine.getDate();

            if (!nextLine.getRunType().equals(runType)) {
              inPostFlushing = (DateTimeUtils.secondsBetween(lineDate, previousDate) <= instrument.getPostFlushingTime());
              finished = true;
            } else if (DateTimeUtils.secondsBetween(lineDate, nextDate) > instrument.getPostFlushingTime()) {
              inPostFlushing = false;
              finished = true;
            } else {
              previousLineIndex = nextLineIndex;
              previousDate = nextDate;
            }
          }
        }
      } catch (DateTimeSpecificationException e) {
        throw new DataFileException(dataFiles.get(fileIndex).getDatabaseId(), line.intValue(), e);
      }
    }

    return inPostFlushing;
  }

  /**
   * Skip the pre-flushing period for the current run type.
   *
   * <p>
   *   After the method is run, {@code firstLine} will point to the first record
   *   after the flushing period, i.e. the first line that can be used. If the Run
   *   Type changes before the pre-flushing time expires, {@code firstLine} will point
   *   to the first line of the new Run Type.
   * </p>
   *
   * <p>
   *   This method assumes that the passed in index is the first line of that
   *   run type.
   * </p>
   *
   * @param fileIndex The file being processed
   * @param runType The run type being processed
   * @param lineIndex The index of the first line in the run type.
   * @throws DataFileException If the file data cannot be processed
   * @return {@code true} if the flushing time expires before the Run Type changes; {@code false} if a new Run Type is encountered
   * @throws FileDefinitionException If an invalid Run Type is found
   */
  private boolean skipPreFlushingTime(int fileIndex, String runType, ExtendedMutableInt lineIndex) throws DataFileException, FileDefinitionException {

    boolean withinRunType = true;

    if (instrument.getPreFlushingTime() > 0) {
      try {
        LocalDateTime lastTime = getLine(fileIndex, lineIndex).getDate();

        lineIndex.increment();
        if (lineIndex.greaterThanOrEqualTo(getFileSize(fileIndex))) {
          lineIndex.setValue(EOF.getValue());
        } else {
          LocalDateTime currentTime = getLine(fileIndex, lineIndex).getDate();
          while (DateTimeUtils.secondsBetween(lastTime, currentTime) < instrument.getPreFlushingTime()) {
            lineIndex.increment();
            if (lineIndex.greaterThanOrEqualTo(getFileSize(fileIndex))) {
              lineIndex.setValue(EOF.getValue());
              break;
            } else {
              DataFileLine newLine = getLine(fileIndex, lineIndex);
              currentTime = newLine.getDate();
              if (!newLine.getRunType().equals(runType)) {
                withinRunType = false;
                break;
              }
            }
          }
        }
      } catch (DateTimeSpecificationException e) {
        throw new DataFileException(dataFiles.get(fileIndex).getDatabaseId(), lineIndex.intValue(), e);
      }

    }

    return withinRunType;
  }

  /**
   * Convenience method to get a line from a file
   * @param fileIndex The file
   * @param line The line index
   * @return The line
   */
  protected DataFileLine getLine(int fileIndex, int line) {
    return data.get(fileIndex).get(line);
  }

  /**
   * Convenience method to get a line from a file
   * @param fileIndex The file
   * @param line The line index
   * @return The line
   */
  protected DataFileLine getLine(int fileIndex, MutableInt line) {
    return getLine(fileIndex, line.getValue());
  }

  /**
   * Convenience method to get number of lines in a file
   * @param fileIndex The file
   * @return The file's size
   */
  protected int getFileSize(int fileIndex) {
    return data.get(fileIndex).size();
  }

  /**
   * Create a new record object from the currently selected lines
   * @return The new record
   * @throws DataFileException If data cannot be extracted
   */
  private DataSetRawDataRecord createRecord() throws DataFileException {

    DataSetRawDataRecord record;

    // TODO This fileId will be set each time a value is read, so any exceptions
    // will contain the correct file ID. This will be written properly for v2.0.0
    long fileId = -1;

    try {
      RunTypeCategory runTypeCategory = getSelectedRunTypeCategory();

      // Only get the position for measurements
      Double longitude = null;
      Double latitude = null;
      if (runTypeCategory.isMeasurementType()) {
       longitude = getSelectedLongitude();
       latitude = getSelectedLatitude();
      }

      record = new DataSetRawDataRecord(dataSet, getSelectedTime(), longitude, latitude, getSelectedRunType(), runTypeCategory);

      for (Map.Entry<SensorType, Set<SensorAssignment>> entry : instrument.getSensorAssignments().entrySet()) {

        SensorType sensorType = entry.getKey();
        Set<SensorAssignment> assignments = entry.getValue();

        // TODO We don't handle diagnostics at this point in the migration.
        // They'll be added back in later
        if (!sensorType.isDiagnostic()) {
          for (SensorAssignment assignment : assignments) {
            Double sensorValue = getSensorValue(assignment);
            record.setSensorValue(sensorType.getName(), sensorValue);
          }
        }

        // This is the old algorithm that does some of the primary/
        // secondary sensor stuff. This may be needed later in the migration.
        /*
        if (instrument.getSensorAssignments().isAssignmentRequired(sensorType)) {

          double primarySensorTotal = 0.0;
          int primarySensorCount = 0;

          double fallbackSensorTotal = 0.0;
          int fallbackSensorCount = 0;

          for (SensorAssignment assignment : assignments) {
            fileId = getFileId(assignment);
            Double sensorValue = getSensorValue(assignment);
            if (null != sensorValue) {
              if (assignment.isPrimary()) {
                primarySensorTotal += sensorValue;
                primarySensorCount++;
              } else {
                fallbackSensorTotal+= sensorValue;
                fallbackSensorCount++;
              }
            }
          }

          Double finalSensorValue = null;

          if (primarySensorCount > 0) {
            finalSensorValue = new Double(primarySensorTotal / primarySensorCount);
          } else if (fallbackSensorCount > 0) {
            finalSensorValue = new Double(fallbackSensorTotal / fallbackSensorCount);
          }

          record.setSensorValue(sensorType.getName(), finalSensorValue);
        } else {
          for (SensorAssignment assignment : assignments) {
            record.setDiagnosticValue(assignment.getDatabaseId(), getSensorValue(assignment));
          }
        }
        */
      }
    } catch (DataFileException e) {
      throw e;
    } catch (Exception e) {
      throw new DataFileException(fileId, DataFileException.NO_LINE_NUMBER, e);
    }

    return record;
  }

  /**
   * Get the date/time for the currently selected line(s)
   * @return The date/time
   * @throws DataFileException If the date/time cannot be extracted
   */
  protected abstract LocalDateTime getSelectedTime() throws DataFileException, DateTimeSpecificationException;

  /**
   * Get the longitude for the currently selected line(s)
   * @return The longitude
   * @throws DataFileException If the file contents cannot be extracted
   * @throws PositionException If the latitude is invalid
   */
  protected abstract double getSelectedLongitude() throws DataFileException, PositionException;

  /**
   * Get the latitude for the currently selected line(s)
   * @return The latitude
   * @throws DataFileException If the file contents cannot be extracted
   * @throws PositionException If the latitude is invalid
   */
  protected abstract double getSelectedLatitude() throws DataFileException, PositionException;

  /**
   * Get the Run Type of the currently selected row(s)
   * @return The Run Type
   * @throws DataFileException If the run type cannot be extracted
   * @throws FileDefinitionException If the run type is invalid
   */
  private String getSelectedRunType() throws DataFileException, FileDefinitionException {
    int fileIndex = getCoreFileIndex();
    int currentRow = selectedRows.get(fileIndex).get(0);
    return data.get(fileIndex).get(currentRow).getRunType();
  }

  /**
   * Get the Run Type Category of the currently selected row(s)
   * @return The Run Type Category
   * @throws DataFileException If the run type cannot be extracted
   * @throws FileDefinitionException If the run type is invalid
   */
  private RunTypeCategory getSelectedRunTypeCategory() throws DataFileException, FileDefinitionException {
    int fileIndex = getCoreFileIndex();
    int currentRow = selectedRows.get(fileIndex).get(0);
    return data.get(fileIndex).get(currentRow).getRunTypeCategory();
  }

  /**
   * Get the primary file definition that will be used as the basis for times
   * and positions in extracted data. This is the first definition that contains
   * Run Types, and therefore the core measurement for the system.
   *
   * @return The index of the first primary file
   */
  protected int getCoreFileIndex() {
    int result = -1;

    for (int i = 0; i < fileDefinitions.size(); i++) {
      if (fileDefinitions.get(i).hasRunTypes()) {
        result = i;
        break;
      }
    }

    return result;
  }

  /**
   * Get the file definition containing position data
   * @return The index of the file containing position data
   */
  protected int getPositionFileIndex() {
    int result = -1;

    for (int i = 0; i < fileDefinitions.size(); i++) {
      if (null != fileDefinitions.get(i).getLatitudeSpecification()) {
        result = i;
        break;
      }
    }

    return result;
  }

  /**
   * Calculate the sensor value for a particular sensor assignment.
   *
   * This will be calculated from all the currently selected lines,
   * with missing values ignored as needed
   *
   * @param assignment The sensor assignment
   * @return The sensor value
   * @throws DataFileException If any field values cannot be extracted
   */
  public Double getSensorValue(SensorAssignment assignment) throws DataFileException {

    int fileIndex = getFileDefinition(assignment.getDataFile());
    List<Integer> rows = selectedRows.get(fileIndex);

    List<Double> rowValues = new ArrayList<Double>();
    for (int row : rows) {
      DataFileLine line = data.get(fileIndex).get(row);
      Double rawValue = line.getFieldValue(assignment.getColumn(), assignment.getMissingValue());
      if (calibrationSet.containsTarget(assignment.getDatabaseId())) {
        rawValue = calibrationSet.getTargetCalibration(String.valueOf(assignment.getDatabaseId()))
            .calibrateValue(rawValue);
      }
      rowValues.add(rawValue);
    }

    return calculateValue(rowValues);
  }

  /**
   * Get the data file ID for a given sensor assignment
   * @param assignment The sensor assignment
   * @return The data file ID
   */
  private long getFileId(SensorAssignment assignment) {
    int fileIndex = getFileDefinition(assignment.getDataFile());
    return dataFiles.get(fileIndex).getDatabaseId();
  }

  /**
   * Calculate the final value from a set of values extracted from
   * the currently selected rows.
   * @param rowValues The row values
   * @return The final calculated value
   */
  protected abstract Double calculateValue(List<Double> rowValues);

  /**
   * Get the index of the file definition with the specified name
   * @param fileName The file definition name
   * @return The file definition index
   */
  private int getFileDefinition(String fileName) {
    int result = -1;

    for (int i = 0; i < fileDefinitions.size(); i++) {
      if (fileDefinitions.get(i).getFileDescription().equals(fileName)) {
        result = i;
        break;
      }
    }

    return result;
  }


  /**
   * Set the calibration set for this data set. The calibration set is used when
   * calculating the calibrated data from the raw data.
   *
   * @param calibrationSet
   */
  private void setCalibrationSet(CalibrationSet calibrationSet) {
    this.calibrationSet = calibrationSet;
  }
}