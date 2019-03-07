package uk.ac.exeter.QuinCe.jobs.files;

import java.sql.Connection;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Properties;
import java.util.Set;

import org.apache.commons.lang3.exception.ExceptionUtils;

import uk.ac.exeter.QuinCe.data.Dataset.DataSet;
import uk.ac.exeter.QuinCe.data.Dataset.DataSetDB;
import uk.ac.exeter.QuinCe.data.Dataset.DataSetDataDB;
import uk.ac.exeter.QuinCe.data.Dataset.InvalidDataSetStatusException;
import uk.ac.exeter.QuinCe.data.Dataset.SensorValue;
import uk.ac.exeter.QuinCe.data.Files.DataFile;
import uk.ac.exeter.QuinCe.data.Files.DataFileDB;
import uk.ac.exeter.QuinCe.data.Instrument.FileDefinition;
import uk.ac.exeter.QuinCe.data.Instrument.Instrument;
import uk.ac.exeter.QuinCe.data.Instrument.InstrumentDB;
import uk.ac.exeter.QuinCe.data.Instrument.SensorDefinition.SensorAssignment;
import uk.ac.exeter.QuinCe.data.Instrument.SensorDefinition.SensorType;
import uk.ac.exeter.QuinCe.jobs.InvalidJobParametersException;
import uk.ac.exeter.QuinCe.jobs.Job;
import uk.ac.exeter.QuinCe.jobs.JobFailedException;
import uk.ac.exeter.QuinCe.jobs.JobThread;
import uk.ac.exeter.QuinCe.utils.DatabaseException;
import uk.ac.exeter.QuinCe.utils.DatabaseUtils;
import uk.ac.exeter.QuinCe.utils.MissingParamException;
import uk.ac.exeter.QuinCe.utils.RecordNotFoundException;
import uk.ac.exeter.QuinCe.web.system.ResourceManager;

/**
 * Job to extract the data for a data set from the uploaded data files
 * @author Steve Jones
 *
 */
public class ExtractDataSetJob extends Job {

  /**
   * The parameter name for the data set id
   */
  public static final String ID_PARAM = "id";

  /**
   * The data set being processed by the job
   */
  private DataSet dataSet = null;

  /**
   * The instrument to which the data set belongs
   */
  private Instrument instrument = null;

  /**
   * Name of the job, used for reporting
   */
  private final String jobName = "Dataset Extraction";

  /**
   * Initialise the job object so it is ready to run
   *
   * @param resourceManager The system resource manager
   * @param config The application configuration
   * @param jobId The id of the job in the database
   * @param parameters The job parameters, containing the file ID
   * @throws InvalidJobParametersException If the parameters are not valid for the job
   * @throws MissingParamException If any of the parameters are invalid
   * @throws RecordNotFoundException If the job record cannot be found in the database
   * @throws DatabaseException If a database error occurs
   */
  public ExtractDataSetJob(ResourceManager resourceManager, Properties config, long jobId, Map<String, String> parameters) throws MissingParamException, InvalidJobParametersException, DatabaseException, RecordNotFoundException {
    super(resourceManager, config, jobId, parameters);
  }

  @Override
  protected void execute(JobThread thread) throws JobFailedException {

    Connection conn = null;

    try {
      conn = dataSource.getConnection();
      conn.setAutoCommit(false);

      // Get the data set from the database
      dataSet = DataSetDB.getDataSet(conn, Long.parseLong(parameters.get(ID_PARAM)));

      // Reset the data set and all associated data
      reset(conn);

      Instrument instrument = InstrumentDB.getInstrument(conn, dataSet.getInstrumentId(),
        ResourceManager.getInstance().getSensorsConfiguration(),
        ResourceManager.getInstance().getRunTypeCategoryConfiguration());

      List<DataFile> files = DataFileDB.getDataFiles(conn,
        ResourceManager.getInstance().getConfig(), dataSet.getSourceFiles(conn));

      List<SensorValue> sensorValues = new ArrayList<SensorValue>();

      for (DataFile file : files) {
        FileDefinition fileDefinition = file.getFileDefinition();

        int currentLine = file.getFirstDataLine();
        while (currentLine < file.getContentLineCount()) {

          List<String> line = file.getLine(currentLine);
          LocalDateTime time = file.getDate(line);

          // Position
          if (null != fileDefinition.getLongitudeSpecification()) {
            double longitude = file.getLongitude(line);
            sensorValues.add(new SensorValue(dataSet.getId(),
              FileDefinition.LONGITUDE_COLUMN_ID, time, String.valueOf(longitude)));
          }

          if (null != fileDefinition.getLongitudeSpecification()) {
            double latitude = file.getLatitude(line);
            sensorValues.add(new SensorValue(dataSet.getId(),
              FileDefinition.LATITUDE_COLUMN_ID, time, String.valueOf(latitude)));
          }

          // Assigned columns
          for (Entry<SensorType, Set<SensorAssignment>> entry :
            instrument.getSensorAssignments().entrySet()) {

            for (SensorAssignment assignment : entry.getValue()) {
              if (assignment.getDataFile().equals(fileDefinition.getFileDescription())) {

                sensorValues.add(new SensorValue(dataSet.getId(),
                  assignment.getDatabaseId(), time,
                  file.getStringValue(line, assignment.getColumn(),
                    assignment.getMissingValue())));
              }
            }
          }

          currentLine++;
        }
      }


      DataSetDataDB.storeSensorValues(conn, sensorValues);

      conn.commit();
    } catch (Exception e) {
      e.printStackTrace();
      DatabaseUtils.rollBack(conn);
      try {
        // Set the dataset to Error status
        dataSet.setStatus(DataSet.STATUS_ERROR);
        // And add a (friendly) message...
        StringBuffer message = new StringBuffer();
        message.append(getJobName());
        message.append(" - error: ");
        message.append(e.getMessage());
        dataSet.addMessage(message.toString(), ExceptionUtils.getStackTrace(e));
        DataSetDB.updateDataSet(conn, dataSet);
        conn.commit();
      } catch (Exception e1) {
        e1.printStackTrace();
      }
      throw new JobFailedException(id, e);
    } finally {
      DatabaseUtils.closeConnection(conn);
    }

    /*
    Connection conn = null;

    try {
      conn = dataSource.getConnection();
      conn.setAutoCommit(false);

      // Get the data set from the database
      dataSet = DataSetDB.getDataSet(conn, Long.parseLong(parameters.get(ID_PARAM)));
      // Reset the data set and all associated data
      reset(conn);

      // Clear messages before executing job
      dataSet.clearMessages();
      // Set processing status
      dataSet.setStatus(DataSet.STATUS_DATA_EXTRACTION);
      DataSetDB.updateDataSet(conn, dataSet);
      conn.commit();

      // Get related data
      instrument = InstrumentDB.getInstrument(conn, dataSet.getInstrumentId(), resourceManager.getSensorsConfiguration(), resourceManager.getRunTypeCategoryConfiguration());

      DataSetRawData rawData = DataSetRawDataFactory.getDataSetRawData(dataSource, dataSet, instrument);

      LocalDateTime realStartTime = null;
      LocalDateTime realEndTime = null;

      DataSetRawDataRecord record = rawData.getNextRecord();
      while (null != record) {

        if (null == realStartTime && null != record.getDate()) {
          realStartTime = record.getDate();
        }

        if (null != record.getDate()) {
          realEndTime = record.getDate();
        }

        if (record.isMeasurement()) {
          DataSetDataDB.storeRecord(conn, record);
        } else if (record.isCalibration()) {
          CalibrationDataDB.storeCalibrationRecord(conn, record);
        }

        // Read the next record
        record = rawData.getNextRecord();
      }

      // Adjust the Dataset limits to the actual extracted data
      if (null != realStartTime) {
        dataSet.setStart(realStartTime);
      }

      if (null != realEndTime) {
        dataSet.setEnd(realEndTime);
      }

      dataSet.setStatus(DataSet.STATUS_DATA_REDUCTION);
      DataSetDB.updateDataSet(conn, dataSet);
      Map<String, String> jobParams = new HashMap<String, String>();
      jobParams.put(DataReductionJob.ID_PARAM, String.valueOf(Long.parseLong(parameters.get(ID_PARAM))));
      JobManager.addJob(dataSource, JobManager.getJobOwner(dataSource, id), DataReductionJob.class.getCanonicalName(), jobParams);

      conn.commit();

    } catch (Exception e) {
      e.printStackTrace();
      DatabaseUtils.rollBack(conn);
      try {
        // Set the dataset to Error status
        dataSet.setStatus(DataSet.STATUS_ERROR);
        // And add a (friendly) message...
        StringBuffer message = new StringBuffer();
        message.append(getJobName());
        message.append(" - error: ");
        message.append(e.getMessage());
        dataSet.addMessage(message.toString(), ExceptionUtils.getStackTrace(e));
        DataSetDB.updateDataSet(conn, dataSet);
        conn.commit();
      } catch (Exception e1) {
        e1.printStackTrace();
      }
      throw new JobFailedException(id, e);
    } finally {
      DatabaseUtils.closeConnection(conn);
    }

    */
  }

  @Override
  protected void validateParameters() throws InvalidJobParametersException {
    // TODO Auto-generated method stub
  }

  /**
   * Reset the data set processing.
   *
   * Delete all related records and reset the status
   *
   * @throws MissingParamException
   *           If any of the parameters are invalid
   * @throws InvalidDataSetStatusException
   *           If the method sets an invalid data set status
   * @throws DatabaseException
   *           If a database error occurs
   * @throws RecordNotFoundException
   *           If the record don't exist
   */
  private void reset(Connection conn)
      throws MissingParamException, InvalidDataSetStatusException,
      DatabaseException, RecordNotFoundException {

    /*
    try {
      CalibrationDataDB.deleteDatasetData(conn, dataSet);
      CalculationDBFactory.getCalculationDB().deleteDatasetCalculationData(conn, dataSet);
      DataSetDB.deleteDatasetData(conn, dataSet);
      dataSet.setStatus(DataSet.STATUS_WAITING);
      DataSetDB.updateDataSet(conn, dataSet);
      conn.commit();
    } catch (SQLException e) {
      throw new DatabaseException("Error while resetting dataset data", e);
    }

    */
  }

  @Override
  public String getJobName() {
    return jobName;
  }
}
