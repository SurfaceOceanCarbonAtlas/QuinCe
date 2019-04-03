package uk.ac.exeter.QuinCe.data.Dataset;

import java.time.LocalDateTime;
import java.util.Collection;

import uk.ac.exeter.QuinCe.data.Dataset.QC.Flag;
import uk.ac.exeter.QuinCe.data.Dataset.QC.Routines.AutoQCResult;
import uk.ac.exeter.QuinCe.data.Dataset.QC.Routines.RoutineException;
import uk.ac.exeter.QuinCe.data.Dataset.QC.Routines.RoutineFlag;
import uk.ac.exeter.QuinCe.utils.RecordNotFoundException;

/**
 * Represents a single sensor value
 * @author Steve Jones
 *
 */
public class SensorValue {

  /**
   * Indicatest that this value has not be stored in the database
   */
  public static final long NO_RECORD = -1;

  /**
   * The database ID of this value
   */
  private final long databaseId;

  /**
   * The ID of the dataset that the sensor value is in
   */
  private final long datasetId;

  /**
   * The ID of the column that the value is in. Either the ID
   * of a row in the {@code file_column} table, or a special value
   * (e.g. for lon/lat)
   */
  private final long columnId;

  /**
   * The time that the value was measured
   */
  private final LocalDateTime time;

  /**
   * The automatic QC result
   */
  private AutoQCResult autoQC = null;

  /**
   * The user QC flag
   */
  private Flag userQCFlag = Flag.ASSUMED_GOOD;

  /**
   * The user QC message
   */
  private String userQCMessage = null;

  /**
   * The value (can be null)
   */
  private final String value;

  /**
   * Indicates whether the value needs to be saved to the database
   */
  private boolean dirty;

  /**
   * Build a sensor value with default QC flags
   * @param datasetId
   * @param columnId
   * @param time
   * @param value
   */
  public SensorValue(long datasetId, long columnId,
    LocalDateTime time, String value) {

    this.databaseId = NO_RECORD;
    this.datasetId = datasetId;
    this.columnId = columnId;
    this.time = time;
    this.value = value;
    this.autoQC = new AutoQCResult();
    this.dirty = true;
  }

  /**
   * Build a sensor value with default QC flags
   * @param datasetId
   * @param columnId
   * @param time
   * @param value
   */
  public SensorValue(long databaseId, long datasetId, long columnId,
    LocalDateTime time, String value, AutoQCResult autoQc,
    Flag userQcFlag, String userQcMessage) {

    this.databaseId = databaseId;
    this.datasetId = datasetId;
    this.columnId = columnId;
    this.time = time;
    this.value = value;

    if (null == autoQc) {
      this.autoQC = new AutoQCResult();
    } else {
      this.autoQC = autoQc;
    }

    this.userQCFlag = userQcFlag;
    this.userQCMessage = userQcMessage;
    this.dirty = false;
  }

  /**
   * Get the database ID of the dataset to which this value belongs
   * @return The dataset ID
   */
  public long getDatasetId() {
    return datasetId;
  }

  /**
   * Get the database ID of the file column from which this value was extracted
   * @return The column ID
   */
  public long getColumnId() {
    return columnId;
  }

  /**
   * Get the time that this value was measured
   * @return The measurement time
   */
  public LocalDateTime getTime() {
    return time;
  }

  /**
   * Get the measured value in its original string format
   * @return The value
   */
  public String getValue() {
    return value;
  }

  /**
   * Get the value as a Double. No error checking is performed.
   * Returns {@code null} if the value is {@code null}.
   *
   * All commas are removed from number before parsing
   *
   * @return The value as a Double
   */
  public Double getDoubleValue() {
    Double result = Double.NaN;
    if (null != value) {
      result = Double.parseDouble(value.replaceAll(",", ""));
    }

    return result;
  }

  /**
   * Indicates whether or not this value is {@code null}.
   * @return {@code true} if the value is null; {@code false} otherwise.
   */
  public boolean isNaN() {
    return getDoubleValue().isNaN();
  }

  /**
   * Get the overall QC flag resulting from the automatic QC.
   * This will be the most significant flag generated by all QC routines
   * run on the value.
   * @return The automatic QC flag
   */
  public Flag getAutoQcFlag() {
    return autoQC.getOverallFlag();
  }

  /**
   * Get the complete automatic QC result
   * @return The automatic QC result
   */
  public AutoQCResult getAutoQcResult() {
    return autoQC;
  }

  /**
   * Get the QC flag set by the user
   * @return The user QC flag
   */
  public Flag getUserQCFlag() {
    return userQCFlag;
  }

  /**
   * Get the QC message entered by the user
   * @return The user QC message
   */
  public String getUserQCMessage() {
    return null == userQCMessage ? "" : userQCMessage;
  }

  /**
   * Reset the automatic QC result
   * @throws RecordNotFoundException If the value has not yet been stored
   *         in the database
   */
  public void clearAutomaticQC() throws RecordNotFoundException {
    if (!isInDatabase()) {
      throw new RecordNotFoundException(
        "SensorValue has not been stored in the database");
    }
    autoQC = new AutoQCResult();

    // Reset the user QC if it hasn't been set by the user
    if (userQCFlag.equals(Flag.ASSUMED_GOOD) || userQCFlag.equals(Flag.NEEDED)) {
      userQCFlag = Flag.ASSUMED_GOOD;
      userQCMessage = null;
    }

    dirty = true;
  }

  /**
   * Add a flag from an automatic QC routine to the automatic QC result
   * @param flag The flag
   * @throws RecordNotFoundException If the value has not yet been stored
   *         in the database
   */
  public void addAutoQCFlag(RoutineFlag flag)
    throws RecordNotFoundException, RoutineException {
    if (!isInDatabase()) {
      throw new RecordNotFoundException(
        "SensorValue has not been stored in the database");
    }
    autoQC.add(flag);

    // Update the user QC if it hasn't been set by the user
    if (userQCFlag.equals(Flag.ASSUMED_GOOD) || userQCFlag.equals(Flag.NEEDED)) {
      userQCFlag = Flag.NEEDED;
      userQCMessage = autoQC.getAllMessages();
    }

    dirty = true;
  }

  /**
   * Set the User QC information
   * @param flag The user QC flag
   * @param message The user QC message
   * @throws RecordNotFoundException If the value has not yet been stored
   *         in the database
   */
  public void setUserQC(Flag flag, String message) throws RecordNotFoundException {
    if (!isInDatabase()) {
      throw new RecordNotFoundException(
        "SensorValue has not been stored in the database");
    }
    userQCFlag = flag;
    userQCMessage = message;
    dirty = true;
  }

  /**
   * Determine whether or not this value needs to be saved to the database
   * @return {@code true} if the value needs to be saved; {@code false} otherwise
   */
  public boolean isDirty() {
    return dirty;
  }

  /**
   * Determine whether or not this value is stored in the database
   * @return {@code true} if the value is in the database; {@code false} if
   *         it is a new record and is yet to be saved
   */
  public boolean isInDatabase() {
    return (databaseId != NO_RECORD);
  }

  /**
   * Clear the {@code dirty} flag on a collection of SensorValues
   * @param sensorValues The values to be cleared
   */
  public static void clearDirtyFlag(Collection<SensorValue> sensorValues) {
    for (SensorValue value : sensorValues) {
      value.dirty = false;
    }
  }

  /**
   * Get the database ID of this sensor value
   * @return The database ID
   */
  public long getDatabaseId() {
    return databaseId;
  }
}