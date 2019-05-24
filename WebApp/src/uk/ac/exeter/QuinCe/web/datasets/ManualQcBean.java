package uk.ac.exeter.QuinCe.web.datasets;

import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import javax.faces.bean.ManagedBean;
import javax.faces.bean.SessionScoped;

import org.primefaces.json.JSONArray;
import org.primefaces.json.JSONObject;

import com.google.gson.Gson;

import uk.ac.exeter.QCRoutines.messages.Flag;
import uk.ac.exeter.QuinCe.data.Calculation.CalculationDBFactory;
import uk.ac.exeter.QuinCe.data.Calculation.CommentSet;
import uk.ac.exeter.QuinCe.data.Calculation.CommentSetEntry;
import uk.ac.exeter.QuinCe.data.Dataset.DataSet;
import uk.ac.exeter.QuinCe.data.Dataset.DataSetDB;
import uk.ac.exeter.QuinCe.data.Dataset.DataSetDataDB;
import uk.ac.exeter.QuinCe.data.Dataset.DiagnosticDataDB;
import uk.ac.exeter.QuinCe.data.Dataset.DataReduction.DataReducerFactory;
import uk.ac.exeter.QuinCe.data.Dataset.QC.InvalidFlagException;
import uk.ac.exeter.QuinCe.data.Dataset.QC.Routines.RoutineException;
import uk.ac.exeter.QuinCe.data.Instrument.FileColumn;
import uk.ac.exeter.QuinCe.data.Instrument.FileDefinition;
import uk.ac.exeter.QuinCe.data.Instrument.InstrumentDB;
import uk.ac.exeter.QuinCe.data.Instrument.SensorDefinition.InstrumentVariable;
import uk.ac.exeter.QuinCe.jobs.JobManager;
import uk.ac.exeter.QuinCe.jobs.files.DataReductionJob;
import uk.ac.exeter.QuinCe.utils.DatabaseException;
import uk.ac.exeter.QuinCe.utils.MissingParamException;
import uk.ac.exeter.QuinCe.web.PlotPageBean;
import uk.ac.exeter.QuinCe.web.Variable;
import uk.ac.exeter.QuinCe.web.VariableList;

/**
 * User QC bean
 * @author Steve Jones
 *
 */
@ManagedBean
@SessionScoped
public class ManualQcBean extends PlotPageBean {

  /**
   * Navigation to the calibration data plot page
   */
  private static final String NAV_PLOT = "user_qc";

  /**
   * Navigation to the dataset list
   */
  private static final String NAV_DATASET_LIST = "dataset_list";

  /**
   * The set of comments for the user QC dialog. Stored as a Javascript array of entries, with each entry containing
   * Comment, Count and Flag value
   */
  private String userCommentList = "[]";

  /**
   * The worst flag set on the selected rows
   */
  private Flag worstSelectedFlag = Flag.GOOD;

  /**
   * The user QC flag
   */
  private int userFlag;

  /**
   * The user QC comment
   */
  private String userComment;

  /**
   * The basic list of table columns. This excludes Date/Time, Lon and Lat,
   * and does not include the sub-columns (i.e. qc info)
   */
  private List<QCColumn> tableColumns;

  /**
   * The table content for the current field set
   */
  private QCTableData tableData;

  /**
   * The list of times in the table data. Used for quick lookups by index
   */
  private List<LocalDateTime> tableTimes;

  /**
   * Initialise the required data for the bean
   */
  @Override
  public void init() {
    setFieldSet(DataSetDataDB.SENSORS_FIELDSET);
  }

  @Override
  protected List<Long> loadRowIds() throws Exception {
    if (null == tableData) {
      loadTableData();
    }

    return Stream.iterate(0L, n -> n + 1)
       .limit(tableData.size() - 1)
       .collect(Collectors.toList());
  }

  @Override
  protected String getScreenNavigation() {
    return NAV_PLOT;
  }

  @Override
  protected String buildPlotLabels(int plotIndex) {
    return null;
  }

  @Override
  protected String loadTableData(int start, int length) throws Exception {

    // TODO Convert to GSON - I don't have the API here

    JSONArray json = new JSONArray();

    for (int i = start; i < start + length; i++) {

      JSONObject obj = new JSONObject();

      obj.put("DT_RowId", tableTimes.get(i));

      int columnIndex = 0;
      obj.put(String.valueOf(columnIndex), tableTimes.get(i));

      LinkedHashMap<Long, QCColumnValue> row = tableData.get(tableTimes.get(i));

      for (QCColumnValue value : row.values()) {

        if (null == value) {
          columnIndex++;
          obj.put(String.valueOf(columnIndex), JSONObject.NULL);

          columnIndex++;
          obj.put(String.valueOf(columnIndex), JSONObject.NULL);

          columnIndex++;
          obj.put(String.valueOf(columnIndex), JSONObject.NULL);

          columnIndex++;
          obj.put(String.valueOf(columnIndex), JSONObject.NULL);

          columnIndex++;
          obj.put(String.valueOf(columnIndex), JSONObject.NULL);
        } else {
          columnIndex++;
          obj.put(String.valueOf(columnIndex), value.getValue());

          columnIndex++;
          obj.put(String.valueOf(columnIndex), value.isUsed());

          columnIndex++;
          obj.put(String.valueOf(columnIndex), value.getQcFlag().getFlagValue());

          columnIndex++;
          obj.put(String.valueOf(columnIndex), value.needsFlag());

          columnIndex++;
          obj.put(String.valueOf(columnIndex), value.getQcComment());
        }
      }

      json.put(obj);
    }

    return json.toString();
  }

  /**
   * Finish the calibration data validation
   * @return Navigation to the data set list
   */
  public String finish() {
    if (dirty) {
      try {
        DataSetDB.setDatasetStatus(getDataSource(), datasetId, DataSet.STATUS_DATA_REDUCTION);
        Map<String, String> jobParams = new HashMap<String, String>();
        jobParams.put(DataReductionJob.ID_PARAM, String.valueOf(datasetId));
        JobManager.addJob(getDataSource(), getUser(), DataReductionJob.class.getCanonicalName(), jobParams);
      } catch (Exception e) {
        e.printStackTrace();
      }
    }
    return NAV_DATASET_LIST;
  }

  /**
   * Apply the automatically generated QC flags to the rows selected in the table
   */
  public void acceptAutoQc() {
    try {
      CalculationDBFactory.getCalculationDB().acceptAutoQc(getDataSource(), getSelectedRowsList());
      dirty = true;
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  /**
   * Returns the set of comments for the WOCE dialog
   * @return The comments for the WOCE dialog
   */
  public String getUserCommentList() {
    return userCommentList;
  }

  /**
   * Dummy: Sets the set of comments for the WOCE dialog
   * @param userCommentList The comments; ignored
   */
  public void setUserCommentList(String userCommentList) {
    // Do nothing
  }

  /**
   * Get the worst flag set on the selected rows
   * @return The worst flag on the selected rows
   */
  public int getWorstSelectedFlag() {
    return worstSelectedFlag.getFlagValue();
  }

  /**
   * Dummy: Set the worst selected flag
   * @param worstSelectedFlag The flag; ignored
   */
  public void setWorstSelectedFlag(int worstSelectedFlag) {
    // Do nothing
  }

  /**
   * Generate the list of comments for the WOCE dialog
   */
  public void generateUserCommentList() {

    worstSelectedFlag = Flag.GOOD;

    JSONArray json = new JSONArray();

    try {
      CommentSet comments = CalculationDBFactory.getCalculationDB().getCommentsForRows(getDataSource(), getSelectedRowsList());
      for (CommentSetEntry entry : comments) {
        json.put(entry.toJson());
        if (entry.getFlag().moreSignificantThan(worstSelectedFlag)) {
          worstSelectedFlag = entry.getFlag();
        }
      }

    } catch (Exception e) {
      e.printStackTrace();
      JSONArray jsonEntry = new JSONArray();
      jsonEntry.put("Existing comments could not be retrieved");
      jsonEntry.put(4);
      jsonEntry.put(1);
      json.put(jsonEntry);
      worstSelectedFlag = Flag.BAD;
    }

    userCommentList = json.toString();
  }

  /**
   * @return the userComment
   */
  public String getUserComment() {
    return userComment;
  }

  /**
   * @param userComment the userComment to set
   */
  public void setUserComment(String userComment) {
    this.userComment = userComment;
  }

  /**
   * @return the userFlag
   */
  public int getUserFlag() {
    return userFlag;
  }

  /**
   * @param userFlag the userFlag to set
   */
  public void setUserFlag(int userFlag) {
    this.userFlag = userFlag;
  }

  /**
   * Apply the entered WOCE flag and comment to the rows selected in the table
   */
  public void applyManualFlag() {
    try {
      CalculationDBFactory.getCalculationDB().applyManualFlag(getDataSource(), getSelectedRowsList(), userFlag, userComment);
      dirty = true;
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  @Override
  protected Variable getDefaultPlot1XAxis() {
    return variables.getVariableWithLabel("Date/Time");
  }

  @Override
  protected List<Variable> getDefaultPlot1YAxis() {
    List<Variable> result = new ArrayList<Variable>(1);
    result.add(variables.getVariableWithLabel("Intake Temperature"));
    return result;
  }

  @Override
  protected Variable getDefaultPlot2XAxis() {
    return variables.getVariableWithLabel("Date/Time");
  }

  @Override
  protected List<Variable> getDefaultPlot2YAxis() {
    List<Variable> result = new ArrayList<Variable>(1);
    result.add(variables.getVariableWithLabel("Final fCO2"));
    return result;
  }

  @Override
  protected Variable getDefaultMap1Variable() {
    return variables.getVariableWithLabel("Intake Temperature");
  }

  @Override
  protected Variable getDefaultMap2Variable() {
    return variables.getVariableWithLabel("Final fCO2");
  }

  @Override
  protected void buildVariableList(VariableList variables) throws Exception {
    DataSetDataDB.populateVariableList(getDataSource(), getDataset(), variables);
    CalculationDBFactory.getCalculationDB().populateVariableList(variables);
    DiagnosticDataDB.populateVariableList(getDataSource(), getDataset(), variables);
  }

  @Override
  protected String getData(List<String> fields) throws Exception {
    return CalculationDBFactory.getCalculationDB().getJsonData(getDataSource(), getDataset(), fields, fields.get(0));
  }

  @Override
  public boolean getHasTwoPlots() {
    return true;
  }

  @Override
  public String getTableHeadings() {
    List<String> headings = new ArrayList<String>(tableColumns.size() * 4 + 12);

    headings.add("Date/Time");

    for (QCColumn column : tableColumns) {
      headings.add(column.getName());
      headings.add(column.getName() + "_used");
      headings.add(column.getName() + "_qc");
      headings.add(column.getName() + "_needsFlag");
      headings.add(column.getName() + "_comment");
    }

    return new Gson().toJson(headings);
  }

  @Override
  public String getFieldSets() {

    TreeMap<Long, List<Integer>> fieldSets = new TreeMap<Long, List<Integer>>();

    // Date/Time
    int col = 0;

    for (QCColumn column : tableColumns) {
      long fieldSet = column.getFieldSet();

      if (null == fieldSets.get(fieldSet)) {
        fieldSets.put(fieldSet, new ArrayList<Integer>());
      }

      fieldSets.get(fieldSet).add(++col);
      fieldSets.get(fieldSet).add(++col);
      fieldSets.get(fieldSet).add(++col);
      fieldSets.get(fieldSet).add(++col);
      fieldSets.get(fieldSet).add(++col);
    }

    return new Gson().toJson(fieldSets);
  }

  /**
   * Load the data for the current field set
   */
  private void loadTableData() throws MissingParamException, DatabaseException,
    InvalidFlagException, RoutineException {

    try {
      tableColumns = new ArrayList<QCColumn>();

      tableColumns.add(new QCColumn(FileDefinition.LONGITUDE_COLUMN_ID,
        "Longitude", 0));
      tableColumns.add(new QCColumn(FileDefinition.LATITUDE_COLUMN_ID,
        "Latitude", 0));

      // Sensor columns
      List<FileColumn> fileColumns = InstrumentDB.getSensorColumns(
        getDataSource(), instrument.getDatabaseId());

      for (FileColumn column : fileColumns) {
        long fieldSet = DataSetDataDB.SENSORS_FIELDSET;
        if (column.getSensorType().isDiagnostic()) {
          fieldSet = DataSetDataDB.DIAGNOSTICS_FIELDSET;
        }

        tableColumns.add(new QCColumn(column.getColumnId(),
          column.getColumnName(), fieldSet));
      }

      // Data reduction columns
      for (InstrumentVariable variable : instrument.getVariables()) {
        LinkedHashMap<String, Long> variableParameters =
          DataReducerFactory.getCalculationParameters(variable);

        // Columns from data reduction are given IDs based on the
        // variable ID and parameter number
        for (Map.Entry<String, Long> entry : variableParameters.entrySet()) {
          tableColumns.add(new QCColumn(entry.getValue(),
            entry.getKey(), variable.getId()));
        }
      }

      tableData = new QCTableData(tableColumns);

      // Load data for sensor columns
      DataSetDataDB.getQCSensorData(getDataSource(), tableData,
        getDataset().getId(), instrument, getTableColumnIDs());

      // Load data reduction data
      DataSetDataDB.getDataReductionData(getDataSource(), tableData, dataset);

      // Extract the times - these will be the row IDs
      tableTimes = new ArrayList<LocalDateTime>(tableData.keySet());
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  private List<Long> getTableColumnIDs() {
    List<Long> result = new ArrayList<Long>(tableColumns.size() + 2);

    for (QCColumn column : tableColumns) {
      result.add(column.getId());
    }

    return result;
  }
}
