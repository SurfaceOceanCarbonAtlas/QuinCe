package uk.ac.exeter.QuinCe.data.Dataset;

import java.time.LocalDateTime;
import java.util.Collection;
import java.util.List;

import uk.ac.exeter.QuinCe.data.Dataset.QC.Flag;
import uk.ac.exeter.QuinCe.data.Dataset.QC.Routines.AutoQCResult;
import uk.ac.exeter.QuinCe.data.Dataset.QC.Routines.RoutineException;
import uk.ac.exeter.QuinCe.data.Dataset.QC.Routines.RoutineFlag;
import uk.ac.exeter.QuinCe.utils.DatabaseUtils;
import uk.ac.exeter.QuinCe.utils.RecordNotFoundException;
import uk.ac.exeter.QuinCe.utils.StringUtils;

/**
 * Represents a single sensor value
 * 
 * @author Steve Jones
 *
 */
public class SensorValue implements Comparable<SensorValue> {

  /**
   * The database ID of this value
   */
  private final long id;

  /**
   * The ID of the dataset that the sensor value is in
   */
  private final long datasetId;

  /**
   * The ID of the column that the value is in. Either the ID of a row in the
   * {@code file_column} table, or a special value (e.g. for lon/lat)
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
   * 
   * @param datasetId
   * @param columnId
   * @param time
   * @param value
   */
  public SensorValue(long datasetId, long columnId, LocalDateTime time,
    String value) {

    this.id = DatabaseUtils.NO_DATABASE_RECORD;
    this.datasetId = datasetId;
    this.columnId = columnId;
    this.time = time;
    this.value = value;
    this.autoQC = new AutoQCResult();
    this.dirty = true;
  }

  /**
   * Build a sensor value with default QC flags
   * 
   * @param datasetId
   * @param columnId
   * @param time
   * @param value
   */
  public SensorValue(long databaseId, long datasetId, long columnId,
    LocalDateTime time, String value, AutoQCResult autoQc, Flag userQcFlag,
    String userQcMessage) {

    this.id = databaseId;
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
   * 
   * @return The dataset ID
   */
  public long getDatasetId() {
    return datasetId;
  }

  /**
   * Get the database ID of the file column from which this value was extracted
   * 
   * @return The column ID
   */
  public long getColumnId() {
    return columnId;
  }

  /**
   * Get the time that this value was measured
   * 
   * @return The measurement time
   */
  public LocalDateTime getTime() {
    return time;
  }

  /**
   * Get the measured value in its original string format
   * 
   * @return The value
   */
  public String getValue() {
    return value;
  }

  /**
   * Get the value as a Double. No error checking is performed. Returns
   * {@code null} if the value is {@code null}.
   *
   * All commas are removed from number before parsing
   *
   * @return The value as a Double
   */
  public Double getDoubleValue() {
    return StringUtils.doubleFromString(value);
  }

  /**
   * Indicates whether or not this value is {@code null}.
   * 
   * @return {@code true} if the value is null; {@code false} otherwise.
   */
  public boolean isNaN() {
    return getDoubleValue().isNaN();
  }

  /**
   * Get the overall QC flag resulting from the automatic QC. This will be the
   * most significant flag generated by all QC routines run on the value.
   * 
   * @return The automatic QC flag
   */
  public Flag getAutoQcFlag() {
    return autoQC.getOverallFlag();
  }

  /**
   * Get the complete automatic QC result
   * 
   * @return The automatic QC result
   */
  public AutoQCResult getAutoQcResult() {
    return autoQC;
  }

  /**
   * Get the QC flag set by the user
   * 
   * @return The user QC flag
   */
  public Flag getUserQCFlag() {
    return userQCFlag;
  }

  /**
   * Get the QC message entered by the user
   * 
   * @return The user QC message
   */
  public String getUserQCMessage() {
    return null == userQCMessage ? "" : userQCMessage;
  }

  /**
   * Reset the automatic QC result
   * 
   * @throws RecordNotFoundException
   *           If the value has not yet been stored in the database
   */
  public void clearAutomaticQC() throws RecordNotFoundException {

    if (!isInDatabase()) {
      throw new RecordNotFoundException(
        "SensorValue has not been stored in the database");
    }
    autoQC = new AutoQCResult();

    // Reset the user QC if it hasn't been set by the user
    if (userQCFlag.equals(Flag.ASSUMED_GOOD)
      || userQCFlag.equals(Flag.NEEDED)) {
      userQCFlag = Flag.ASSUMED_GOOD;
      userQCMessage = null;
    }

    dirty = true;
  }

  /**
   * Add a flag from an automatic QC routine to the automatic QC result
   * 
   * @param flag
   *          The flag
   * @throws RecordNotFoundException
   *           If the value has not yet been stored in the database
   */
  public void addAutoQCFlag(RoutineFlag flag)
    throws RecordNotFoundException, RoutineException {

    if (!isInDatabase()) {
      throw new RecordNotFoundException(
        "SensorValue has not been stored in the database");
    }
    autoQC.add(flag);

    // Update the user QC if it hasn't been set by the user
    if (userQCFlag.equals(Flag.ASSUMED_GOOD)
      || userQCFlag.equals(Flag.NEEDED)) {
      userQCFlag = Flag.NEEDED;
      userQCMessage = autoQC.getAllMessages();
    }

    dirty = true;
  }

  /**
   * Set the User QC information
   * 
   * @param flag
   *          The user QC flag
   * @param message
   *          The user QC message
   * @throws RecordNotFoundException
   *           If the value has not yet been stored in the database
   */
  public void setUserQC(Flag flag, String message)
    throws RecordNotFoundException {
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
   * 
   * @return {@code true} if the value needs to be saved; {@code false}
   *         otherwise
   */
  public boolean isDirty() {
    return dirty;
  }

  /**
   * Determine whether or not this value is stored in the database
   * 
   * @return {@code true} if the value is in the database; {@code false} if it
   *         is a new record and is yet to be saved
   */
  public boolean isInDatabase() {
    return (id != DatabaseUtils.NO_DATABASE_RECORD);
  }

  /**
   * Clear the {@code dirty} flag on a collection of SensorValues
   * 
   * @param sensorValues
   *          The values to be cleared
   */
  public static void clearDirtyFlag(Collection<SensorValue> sensorValues) {
    for (SensorValue value : sensorValues) {
      value.dirty = false;
    }
  }

  /**
   * Get the database ID of this sensor value
   * 
   * @return The database ID
   */
  public long getId() {
    return id;
  }

  @Override
  public int compareTo(SensorValue o) {
    // If the IDs are equal, the objects are equal.
    // Otherwise compare on time, dataset ID, column ID

    int result = 0;

    if (id > -1) {
      result = Long.compare(id, o.id);
    }

    if (result == 0) {
      result = time.compareTo(o.time);
    }

    if (result == 0) {
      result = Long.compare(datasetId, o.datasetId);
    }

    if (result == 0) {
      result = Long.compare(columnId, o.columnId);
    }

    return result;
  }

  /**
   * Clear the automatic QC information for a set of SensorValues
   * 
   * @param values
   *          The values
   * @throws RecordNotFoundException
   */
  public static void clearAutoQC(List<SensorValue> values)
    throws RecordNotFoundException {
    for (SensorValue value : values) {
      value.clearAutomaticQC();
    }
  }

  public static boolean contains(List<SensorValue> values, String searchValue) {
    boolean result = false;

    for (SensorValue testValue : values) {
      if (testValue.getValue().equals(searchValue)) {
        result = true;
        break;
      }
    }

    return result;
  }

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + (int) (columnId ^ (columnId >>> 32));
    result = prime * result + (int) (datasetId ^ (datasetId >>> 32));
    result = prime * result + ((time == null) ? 0 : time.hashCode());
    return result;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj)
      return true;
    if (obj == null)
      return false;
    if (!(obj instanceof SensorValue))
      return false;
    SensorValue other = (SensorValue) obj;
    if (columnId != other.columnId)
      return false;
    if (datasetId != other.datasetId)
      return false;
    if (time == null) {
      if (other.time != null)
        return false;
    } else if (!time.equals(other.time))
      return false;
    return true;
  }
}
