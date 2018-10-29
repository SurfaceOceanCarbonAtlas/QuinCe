package uk.ac.exeter.QuinCe.data.Dataset;

import java.sql.Connection;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

import org.primefaces.json.JSONArray;

import uk.ac.exeter.QuinCe.data.Files.DataFileDB;
import uk.ac.exeter.QuinCe.utils.DatabaseException;
import uk.ac.exeter.QuinCe.utils.DatabaseUtils;
import uk.ac.exeter.QuinCe.utils.DateTimeUtils;
import uk.ac.exeter.QuinCe.utils.Message;
import uk.ac.exeter.QuinCe.utils.MissingParamException;

/**
 * Object to represent a data set
 * @author Steve Jones
 *
 */
public class DataSet {

  /**
   * The numeric value for the error status.
   * The data set will be given this status whenever a processing job fails.
   */
  public static final int STATUS_ERROR = -1;

  /**
   * The string for the error status
   */
  public static final String STATUS_ERROR_NAME = "ERROR";

  /**
   * The numeric value for the data extraction status
   */
  public static final int STATUS_WAITING = 0;

  /**
   * The string for the data extraction status
   */
  public static final String STATUS_WAITING_NAME = "Waiting";

  /**
   * The numeric value for the data extraction status
   */
  public static final int STATUS_DATA_EXTRACTION = 1;

  /**
   * The string for the data extraction status
   */
  public static final String STATUS_DATA_EXTRACTION_NAME = "Data extraction";

  /**
   * The numeric value for the data reduction status
   */
  public static final int STATUS_DATA_REDUCTION = 2;

  /**
   * The string for the data extraction status
   */
  public static final String STATUS_DATA_REDUCTION_NAME = "Data reduction";

  /**
   * The numeric value for the automatic QC status
   */
  public static final int STATUS_AUTO_QC = 3;

  /**
   * The string for the automatic QC status
   */
  public static final String STATUS_AUTO_QC_NAME = "Automatic QC";

  /**
   * The numeric value for the user QC status
   */
  public static final int STATUS_USER_QC = 4;

  /**
   * The string for the user QC status
   */
  public static final String STATUS_USER_QC_NAME = "Ready for QC";

  /**
   * The numeric value for the ready for submission status
   */
  public static final int STATUS_READY_FOR_SUBMISSION = 5;

  /**
   * The string for the ready for submission status
   */
  public static final String STATUS_READY_FOR_SUBMISSION_NAME = "Ready for Submission";

  /**
   * The numeric value for the waiting for approval status
   */
  public static final int STATUS_WAITING_FOR_APPROVAL = 6;

  /**
   * The string for the waiting for approval status
   */
  public static final String STATUS_WAITING_FOR_APPROVAL_NAME = "Waiting for Approval";

  /**
   * The numeric value for the ready for export status
   */
  public static final int STATUS_READY_FOR_EXPORT = 7;

  /**
   * The string for the ready for export status
   */
  public static final String STATUS_READY_FOR_EXPORT_NAME = "Waiting for Automatic Export";

  /**
   * The numeric value for the exporting status
   */
  public static final int STATUS_EXPORTING = 8;

  /**
   * The string for the exporting status
   */
  public static final String STATUS_EXPORTING_NAME = "Automatic export in progress";

  /**
   * The numeric value for the export complete status
   */
  public static final int STATUS_EXPORT_COMPLETE = 9;

  /**
   * The string for the export complete status
   */
  public static final String STATUS_EXPORT_COMPLETE_NAME = "Automatic export complete";

  /**
   * The database ID
   */
  private long id = DatabaseUtils.NO_DATABASE_RECORD;

  /**
   * The database ID of the instrument to which this data set belongs
   */
  private long instrumentId;

  /**
   * The data set name
   */
  private String name;

  /**
   * The start date of the data set
   */
  private LocalDateTime start;

  /**
   * The end date of the data set
   */
  private LocalDateTime end;

  /**
   * Additional properties of the data set
   */
  private Properties properties;

  /**
   * The time when the data file was last touched
   */
  private LocalDateTime lastTouched;

  /**
   * The data set's status
   */
  private int status = STATUS_WAITING;

  /**
   * The date that the status was set
   */
  private LocalDateTime statusDate = null;

  /**
   * Indicates whether or not this is a NRT dataset
   */
  private boolean nrt = false;

  /**
   * Messages from jobs handling this data set
   */
  private ArrayList<Message> messages = new ArrayList<Message>();

  /**
   * The number of records that need flagging
   */
  private int needsFlagCount = -1;

  /**
   * Constructor for all fields
   * @param id
   *    Data set's database ID
   * @param instrumentId
   *    Database ID of the instrument to which the data set belongs
   * @param name
   *    Dataset name
   * @param start
   *    Start date
   * @param end
   *    End date
   * @param status
   *          The current status
   * @param status
   *          The date that the status was set
   * @param nrt
   *    Indicates whether or not this is a NRT dataset
   * @param properties
   *    Additional properties
   * @param lastTouched
   *    Date that the dataset was last accessed
   * @param needsFlagCount
   *    Number of records that need flagging by the user
   * @param messages
   *    List of messages concerning the dataset (errors etc)
   */
  protected DataSet(long id, long instrumentId, String name,
      LocalDateTime start, LocalDateTime end, int status, LocalDateTime statusDate,
      boolean nrt, Properties properties, LocalDateTime lastTouched, int needsFlagCount,
      List<Message> messages) {
    this.id = id;
    this.instrumentId = instrumentId;
    this.name = name;
    this.start = start;
    this.end = end;
    this.status = status;
    this.statusDate = statusDate;
    this.nrt = nrt;
    this.properties = properties;
    this.lastTouched = lastTouched;
    this.needsFlagCount = needsFlagCount;
    this.messages = new ArrayList<Message>(messages);
  }

  /**
   * Constructor for a new, empty data set
   * @param instrumentId The database ID of the instrument to which the data set belongs
   */
  public DataSet(long instrumentId) {
    this.instrumentId = instrumentId;
    this.statusDate = DateTimeUtils.longToDate(System.currentTimeMillis());
  }

  /**
   * Constructor for all fields
   * @param id
   *    Data set's database ID
   * @param instrumentId
   *    Database ID of the instrument to which the data set belongs
   * @param name
   *    Dataset name
   * @param start
   *    Start date
   * @param end
   *    End date
   * @param nrt
   *    Indicates whether or not this is a NRT dataset
   */
  public DataSet(long instrumentId, String name,
      LocalDateTime start, LocalDateTime end, boolean nrt) {
    this.instrumentId = instrumentId;
    this.name = name;
    this.start = start;
    this.end = end;
    this.nrt = nrt;
    this.statusDate = DateTimeUtils.longToDate(System.currentTimeMillis());
  }

  /**
   * Get the data set's status
   * @return The status
   */
  public int getStatus() {
    return status;
  }

  /**
   * Get the human-readable status of the data set
   * @return The status
   */
  public String getStatusName() {
    return getStatusName(status);
  }

  /**
   * Get the human-readable name of a given data set status
   * @param statusValue The status value
   * @return The status name
   */
  public static String getStatusName(int statusValue) {
    String result;

    switch (statusValue) {
    case STATUS_ERROR: {
      result = STATUS_ERROR_NAME;
      break;
    }
    case STATUS_WAITING: {
      result = STATUS_WAITING_NAME;
      break;
    }
    case STATUS_DATA_EXTRACTION: {
      result = STATUS_DATA_EXTRACTION_NAME;
      break;
    }
    case STATUS_DATA_REDUCTION: {
      result = STATUS_DATA_REDUCTION_NAME;
      break;
    }
    case STATUS_AUTO_QC: {
      result = STATUS_AUTO_QC_NAME;
      break;
    }
    case STATUS_USER_QC: {
      result = STATUS_USER_QC_NAME;
      break;
    }
    case STATUS_READY_FOR_SUBMISSION: {
      result = STATUS_READY_FOR_SUBMISSION_NAME;
      break;
    }
    case STATUS_WAITING_FOR_APPROVAL: {
      result = STATUS_WAITING_FOR_APPROVAL_NAME;
      break;
    }
    case STATUS_READY_FOR_EXPORT: {
      result = STATUS_READY_FOR_EXPORT_NAME;
      break;
    }
    case STATUS_EXPORTING: {
      result = STATUS_EXPORTING_NAME;
      break;
    }
    case STATUS_EXPORT_COMPLETE: {
      result = STATUS_EXPORT_COMPLETE_NAME;
      break;
    }
    default: {
      result = "UNKNOWN";
    }
    }

    return result;
  }

  /**
   * Get the data set's database ID
   * @return The database ID
   */
  public long getId() {
    return id;
  }

  /**
   * Set the data set's database ID
   * @param id The database ID
   */
  protected void setId(long id) {
    this.id = id;
  }

  /**
   * Get the database ID of the instrument to which this data set belongs
   * @return The instrument's database ID
   */
  public long getInstrumentId() {
    return instrumentId;
  }

  /**
   * Get the name of the data set
   * @return The data set name
   */
  public String getName() {
    return name;
  }

  /**
   * Set the name of the data set
   * @param name The data set name
   */
  public void setName(String name) {
    this.name = name;
  }

  /**
   * Get the start date of the data set
   * @return The start date
   */
  public LocalDateTime getStart() {
    return start;
  }

  /**
   * Set the start date of the data set
   * @param start The start date
   */
  public void setStart(LocalDateTime start) {
    this.start = start;
  }

  /**
   * Get the end date of the data set
   * @return The end date
   */
  public LocalDateTime getEnd() {
    return end;
  }

  /**
   * Set the end date of the data set
   * @param end The end date
   */
  public void setEnd(LocalDateTime end) {
    this.end = end;
  }

  /**
   * Get the date that the data set was last accessed
   * @return The last access date
   */
  public LocalDateTime getLastTouched() {
    return lastTouched;
  }

  /**
   * Set the last access date of the data set to the current time
   */
  public void touch() {
    this.lastTouched = DateTimeUtils.longToDate(System.currentTimeMillis());
  }

  /**
   * Set the data set's status
   * @param status The status
   * @throws InvalidDataSetStatusException If the status is invalid
   */
  public void setStatus(int status) throws InvalidDataSetStatusException {
    if (!validateStatus(status)) {
      throw new InvalidDataSetStatusException(status);
    }

    this.status = status;
    statusDate = DateTimeUtils.longToDate(System.currentTimeMillis());
  }

  /**
   * Get the date that the dataset's status was set
   * @return The status date
   */
  public LocalDateTime getStatusDate() {
    return statusDate;
  }

  /**
   * Set a property on the data set.
   * @param key The key
   * @param value The value
   * @see Properties#setProperty(String, String)
   */
  public void setProperty(String key, String value) {
    if (null == properties) {
      properties = new Properties();
    }

    properties.setProperty(key, value);
  }

  /**
   * Get a property from the data set
   * @param key The key
   * @return The value
   * @see Properties#getProperty(String)
   */
  public String getProperty(String key) {
    String result = null;

    if (null != properties) {
      result = properties.getProperty(key);
    }

    return result;
  }

  /**
   * Determine whether or not a given status value is valid
   * @param status The status to be checked
   * @return {@code true} if the status is valid; {@code false} if it is not
   */
  public static boolean validateStatus(int status) {
    return (status >= STATUS_ERROR && status <= STATUS_EXPORT_COMPLETE);
  }

  /**
   * Get the number of records that need flagging by the user.
   *
   * This only gives a value if the dataset's status is
   * {@line #STATUS_USER_QC} or later;
   * otherwise it returns {@code -1}.
   * @return The number of records that need flagging
   */
  public int getNeedsFlagCount() {
    int result = needsFlagCount;

    // Status values increase through processing.
    // This value isn't available before automatic
    // QC is finished
    if (getStatus() < STATUS_USER_QC) {
      result = -1;
    }

    return result;
  }

  /**
   * Determine whether or not this dataset can be exported.
   *
   * <p>The data set can be exported if it meets the following
   * criteria:</p>
   * <ul>
   *   <li>The dataset ready for user QC (not in any other state)</li>
   *   <li>No rows need flagging</li>
   * </ul>
   *
   * @return {@code true} if the dataset can be exported; {@code false} if it cannot
   */
  public boolean getCanBeExported() {
    return (getStatus() >= STATUS_USER_QC && needsFlagCount == 0);
  }

  public void addMessage(String message, String details) {
    messages.add(new Message(message, details));
  }

  public List<Message> getMessages() {
    return messages;
  }

  public String getMessagesAsJSONString() {
    JSONArray json = new JSONArray();
    for (Message message : getMessages()) {
      json.put(message.getAsJSON());
    }
    return json.toString();
  }

  public int getMessageCount() {
    return messages.size();
  }

  /**
   * Remove all messages from this dataset.
   */
  public void clearMessages() {
    messages.clear();
  }

  /**
   * Get a list of the raw data files used to construct
   * this DataSet
   * @param dataSource A data source
   * @return The IDs of the files
   * @throws MissingParamException If any required parameters are missing
   * @throws DatabaseException If a database error occurs
   */
  public List<Long> getSourceFiles(Connection conn) throws MissingParamException, DatabaseException {
    return DataFileDB.getFilesWithinDates(conn, instrumentId, start, end);
  }

  /**
   * Determine whether or not this is a NRT dataset
   * @return {@code true} if this is an NRT dataset; {@code false} if it is not
   */
  public boolean isNrt() {
    return nrt;
  }

  @Override
  public boolean equals(Object o) {
    // Compares database ID only
    boolean result = false;

    if (o instanceof DataSet) {
      result = ((DataSet) o).id == id;
    }

    return result;
  }
}
