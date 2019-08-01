package uk.ac.exeter.QuinCe.web.datasets.export;

import java.io.ByteArrayOutputStream;
import java.io.OutputStream;
import java.math.RoundingMode;
import java.sql.Connection;
import java.text.DecimalFormat;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Properties;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import javax.faces.bean.ManagedBean;
import javax.faces.bean.SessionScoped;
import javax.faces.context.ExternalContext;
import javax.faces.context.FacesContext;
import javax.sql.DataSource;

import org.primefaces.json.JSONArray;
import org.primefaces.json.JSONObject;

import uk.ac.exeter.QCRoutines.config.InvalidDataTypeException;
import uk.ac.exeter.QCRoutines.data.NoSuchColumnException;
import uk.ac.exeter.QCRoutines.messages.MessageException;
import uk.ac.exeter.QuinCe.data.Dataset.DataSet;
import uk.ac.exeter.QuinCe.data.Dataset.DataSetDB;
import uk.ac.exeter.QuinCe.data.Dataset.DataSetDataDB;
import uk.ac.exeter.QuinCe.data.Dataset.DataReduction.DataReductionException;
import uk.ac.exeter.QuinCe.data.Export.ExportConfig;
import uk.ac.exeter.QuinCe.data.Export.ExportException;
import uk.ac.exeter.QuinCe.data.Export.ExportOption;
import uk.ac.exeter.QuinCe.data.Files.DataFile;
import uk.ac.exeter.QuinCe.data.Files.DataFileDB;
import uk.ac.exeter.QuinCe.data.Instrument.FileDefinition;
import uk.ac.exeter.QuinCe.data.Instrument.Instrument;
import uk.ac.exeter.QuinCe.utils.DatabaseException;
import uk.ac.exeter.QuinCe.utils.DatabaseUtils;
import uk.ac.exeter.QuinCe.utils.DateTimeUtils;
import uk.ac.exeter.QuinCe.utils.MissingParamException;
import uk.ac.exeter.QuinCe.utils.RecordNotFoundException;
import uk.ac.exeter.QuinCe.utils.StringUtils;
import uk.ac.exeter.QuinCe.web.BaseManagedBean;
import uk.ac.exeter.QuinCe.web.datasets.data.Field;
import uk.ac.exeter.QuinCe.web.datasets.data.FieldValue;
import uk.ac.exeter.QuinCe.web.system.ResourceManager;

@ManagedBean
@SessionScoped
public class ExportBean extends BaseManagedBean {

  /**
   * Navigation to the export page
   */
  private static final String NAV_EXPORT_PAGE = "export";

  /**
   * The database ID of the dataset to be exported
   */
  private DataSet dataset = null;

  /**
   * The chosen export option
   */
  private int chosenExportOption = -1;

  /**
   * Indicates whether raw files should be included in the export
   */
  private boolean includeRawFiles = false;

  /**
   * Formatter for numeric values
   * All values are displayed to 3 decimal places.
   */
  private static DecimalFormat numberFormatter;

  static {
    numberFormatter = new DecimalFormat("#0.000");
    numberFormatter.setRoundingMode(RoundingMode.HALF_UP);
  }

  /**
   * Initialise the bean
   */
  public String start() {
    return NAV_EXPORT_PAGE;
  }

  /**
   * Get the dataset ID
   * @return The dataset ID
   */
  public long getDatasetId() {
    long result = -1;
    if (dataset != null) {
      result = dataset.getId();
    }

    return result;
  }

  /**
   * Set the dataset using its ID
   * @param datasetId The dataset ID
   */
  public void setDatasetId(long datasetId) {
    try {
      this.dataset = DataSetDB.getDataSet(getDataSource(), datasetId);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  /**
   * Get the dataset
   * @return The dataset
   */
  public DataSet getDataset() {
    return dataset;
  }

  /**
   * Set the dataset
   * @param dataset The dataset
   */
  public void setDataset(DataSet dataset) {
    this.dataset = dataset;
  }

  /**
   * Get the list of available file export options
   * @return The export options
   * @throws ExportException In an error occurs while retrieving the export options
   */
  public List<ExportOption> getExportOptions() throws ExportException {
    List<ExportOption> options = ExportConfig.getInstance().getOptions();
    if (chosenExportOption == -1 && options.size() > 0) {
      chosenExportOption = options.get(0).getIndex();
    }

    return options;
  }

  /**
   * Return the ID of the chosen file export option
   * @return The export option ID
   */
  public int getChosenExportOption() {
    return chosenExportOption;
  }

  /**
   * Set the ID of the chosen export option
   * @param chosenExportOption The export option ID
   */
  public void setChosenExportOption(int chosenExportOption) {
    this.chosenExportOption = chosenExportOption;
  }

  /**
   * Export the dataset in the chosen format
   */
  public void exportDataset() {

    if (includeRawFiles) {
      exportDatasetWithRawFiles();
    } else {
      exportSingleFile();
    }
  }

  /**
   * Export a dataset file by itself
   */
  private void exportSingleFile() {

    Connection conn = null;

    try {
      conn = getDataSource().getConnection();
      ExportOption exportOption = getExportOptions().get(chosenExportOption);

      byte[] fileContent = getDatasetExport(conn, getCurrentInstrument(), dataset, exportOption);

      FacesContext fc = FacesContext.getCurrentInstance();
      ExternalContext ec = fc.getExternalContext();

      ec.responseReset();
      ec.setResponseContentType("text/csv");

      // Set it with the file size. This header is optional. It will work if it's omitted,
      // but the download progress will be unknown.
      ec.setResponseContentLength(fileContent.length);

      // The Save As popup magic is done here. You can give it any file name you want, this only won't work in MSIE,
      // it will use current request URL as file name instead.
      ec.setResponseHeader("Content-Disposition", "attachment; filename=\"" + getExportFilename(exportOption) + "\"");

      OutputStream outputStream = ec.getResponseOutputStream();
      outputStream.write(fileContent);

      fc.responseComplete();
    } catch (Exception e) {
      e.printStackTrace();
    } finally {
      DatabaseUtils.closeConnection(conn);
    }
  }

  /**
   * Create a ZIP file containing a dataset and its source raw data files
   */
  private void exportDatasetWithRawFiles() {

    Connection conn = null;

    try {
      conn = getDataSource().getConnection();
      ExportOption exportOption = getExportOptions().get(chosenExportOption);

      byte[] outBytes = buildExportZip(conn, getCurrentInstrument(), dataset, exportOption);

      FacesContext fc = FacesContext.getCurrentInstance();
      ExternalContext ec = fc.getExternalContext();

      ec.responseReset();
      ec.setResponseContentType("application/zip");

      // Set it with the file size. This header is optional. It will work if it's omitted,
      // but the download progress will be unknown.
      ec.setResponseContentLength(outBytes.length);

      // The Save As popup magic is done here. You can give it any file name you want, this only won't work in MSIE,
      // it will use current request URL as file name instead.
      ec.setResponseHeader("Content-Disposition", "attachment; filename=\"" + dataset.getName() + ".zip\"");

      OutputStream outputStream = ec.getResponseOutputStream();
      outputStream.write(outBytes);

      fc.responseComplete();
    } catch (Exception e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } finally {
      DatabaseUtils.closeConnection(conn);
    }
  }

  /**
   * Obtain a dataset in the specified format as a byte array ready for export or
   * storage
   * @param dataSource A data source
   * @param dataSet The dataset
   * @param exportOption The export format
   * @return The exported dataset
   * @throws NoSuchColumnException
   * @throws InvalidDataTypeException
   * @throws MissingParamException
   * @throws DatabaseException
   * @throws RecordNotFoundException
   * @throws MessageException
   * @throws DataReductionException
   */
  private static byte[] getDatasetExport(Connection conn, Instrument instrument, DataSet dataset, ExportOption exportOption)
      throws Exception {


    DataSource dataSource = ResourceManager.getInstance().getDBDataSource();

    ExportData data = new ExportData(
      dataSource, instrument, exportOption);

    // Load sensor data
    List<Long> fieldIds = new ArrayList<Long>();
    fieldIds.add(FileDefinition.LONGITUDE_COLUMN_ID);
    fieldIds.add(FileDefinition.LATITUDE_COLUMN_ID);
    fieldIds.addAll(instrument.getSensorAssignments().getFileColumnIDs());

    DataSetDataDB.getQCSensorData(dataSource, data,
      dataset.getId(), instrument, fieldIds);

    // Load data reduction data
    DataSetDataDB.getDataReductionData(dataSource, data, dataset);


    // Let's make some output
    StringBuilder output = new StringBuilder();

    // Headers
    List<String> headers = new ArrayList<String>();
    headers.add("Date/Time");
    headers.addAll(data.getFieldSets().getTableHeadingsList());
    output.append(StringUtils.collectionToDelimited(headers, exportOption.getSeparator()));
    output.append('\n');

    for (LocalDateTime rowId : data.keySet()) {
      LinkedHashMap<Field, FieldValue> row = data.get(rowId);

      output.append(DateTimeUtils.toIsoDate(rowId));

      Iterator<FieldValue> valueIter = row.values().iterator();
      while (valueIter.hasNext()) {
        output.append(exportOption.getSeparator());
        FieldValue fieldValue = valueIter.next();
        if (null == fieldValue || fieldValue.isNaN()) {
          output.append(exportOption.getMissingValue());
        } else {
          output.append(numberFormatter.format(fieldValue.getValue()));
        }
      }

      output.append('\n');
    }

/*


    // TODO This will get all sensor columns. When the sensor data storage is updated (Issue #576), this can be revised.
    List<DataSetRawDataRecord> datasetData = DataSetDataDB.getMeasurements(conn, dataset);

    List<CalculationRecord> calculationData = new ArrayList<CalculationRecord>(datasetData.size());
    for (DataSetRawDataRecord record : datasetData) {
      CalculationRecord calcRecord = CalculationRecordFactory.makeCalculationRecord(dataset.getId(), record.getId());
      CalculationDBFactory.getCalculationDB().loadCalculationValues(conn, calcRecord);
      calculationData.add(calcRecord);
    }

    StringBuilder output = new StringBuilder();

    // The header
    output.append("Timestamp");
    output.append(exportOption.getSeparator());
    output.append("Longitude");
    output.append(exportOption.getSeparator());
    output.append("Latitude");
    output.append(exportOption.getSeparator());

    for (String sensorColumnHeading : exportOption.getSensorColumnHeadings()) {
      output.append(sensorColumnHeading);
      output.append(exportOption.getSeparator());
    }

    // TODO Replace when mutiple calculation paths are in place
    List<String> calculationColumnHeadings = exportOption.getCalculationColumnHeadings("equilibrator_pco2");
    if (null != calculationColumnHeadings) {
      for (int i = 0; i < calculationColumnHeadings.size(); i++) {
        output.append(calculationColumnHeadings.get(i));
        output.append(exportOption.getSeparator());
      }
    }

    output.append("QC Flag");
    output.append(exportOption.getSeparator());
    output.append("QC Message");
    output.append('\n');

    for (int i = 0; i < datasetData.size(); i++) {

      DataSetRawDataRecord sensorRecord = datasetData.get(i);
      CalculationRecord calculationRecord = calculationData.get(i);

      Flag recordFlag;
      if (dataset.isNrt()) {
        recordFlag = calculationRecord.getAutoFlag();
      } else {
        recordFlag = calculationRecord.getUserFlag();
      }

      if (exportOption.flagAllowed(recordFlag)) {

        output.append(DateTimeUtils.toIsoDate(sensorRecord.getDate()));
        output.append(exportOption.getSeparator());
        output.append(numberFormatter.format(sensorRecord.getLongitude()));
        output.append(exportOption.getSeparator());
        output.append(numberFormatter.format(sensorRecord.getLatitude()));
        output.append(exportOption.getSeparator());

        for (String sensorColumn : exportOption.getSensorColumns()) {

          Double value = null;

          // This is a horrible hack. Equilibrator Pressure should be dealt with as a calculated
          // value. This will be sorted out by issue #1049.
          if (sensorColumn.equals("Equilibrator Pressure")) {
            Double absoluteEquilibratorPressure = sensorRecord.getSensorValue("Equilibrator Pressure (absolute)");
            if (null != absoluteEquilibratorPressure) {
              value = absoluteEquilibratorPressure;
            } else {
              Double differential = sensorRecord.getSensorValue("Equilibrator Pressure (differential)");
              Double atmospheric = sensorRecord.getSensorValue("Ambient Pressure");
              if (null != differential && null != atmospheric) {
                value = atmospheric + differential;
              }
            }
          } else {
            value = sensorRecord.getSensorValue(sensorColumn);
          }

          if (null == value) {
            output.append("NaN");
          } else {
            output.append(numberFormatter.format(value));
          }

          output.append(exportOption.getSeparator());
        }

        List<String> calculationColumns = exportOption.getCalculationColumns("equilibrator_pco2");
        if (null != calculationColumns) {
          for (String calculatedColumn : calculationColumns) {
            Double value = calculationRecord.getNumericValue(calculatedColumn);
            if (null == value) {
              output.append(MISSING_VALUE);
            } else {
              output.append(numberFormatter.format(value));
            }

            output.append(exportOption.getSeparator());
          }
        }

        if (dataset.isNrt()) {
          output.append(Flag.getWoceValue(calculationRecord.getAutoFlag().getFlagValue()));
        } else {
          output.append(Flag.getWoceValue(calculationRecord.getUserFlag().getFlagValue()));
        }

        output.append(exportOption.getSeparator());

        String qcMessage = null;
        if (dataset.isNrt()) {
          qcMessage = calculationRecord.getAutoQCMessagesString();
        } else {
          qcMessage = calculationRecord.getUserMessage();
        }

        if (null != qcMessage) {
          if (qcMessage.length() > 0) {
            output.append(StringUtils.makeCsvString(qcMessage.trim()));
          }
        }

        output.append('\n');
      }

    }
*/

    return output.toString().getBytes();
  }

  /**
   * Get the filename of the file that will be exported
   * @param exportOption The export option
   * @return The export filename
   * @throws Exception If any errors occur
   */
  private String getExportFilename(ExportOption exportOption) throws Exception {
    StringBuffer fileName = new StringBuffer(dataset.getName().replaceAll("\\.", "_"));
    fileName.append('-');
    fileName.append(exportOption.getName());

    if (exportOption.getSeparator().equals("\t")) {
      fileName.append(".tsv");
    } else {
      fileName.append(".csv");
    }

    return fileName.toString();
  }

  /**
   * Determine whether raw files should be included in the export
   * @return {@code true} if raw files should be included; {@code false} if not
   */
  public boolean getIncludeRawFiles() {
    return includeRawFiles;
  }

  /**
   * Specify whether raw files should be included in the export
   * @param includeRawFiles The raw files flag
   */
  public void setIncludeRawFiles(boolean includeRawFiles) {
    this.includeRawFiles = includeRawFiles;
  }

  /**
   * Create the contents of the manifest.json file
   * @return The manifest JSON
   */
  private static JSONObject makeManifest(Connection conn, DataSet dataset,
      List<ExportOption> exportOptions, List<DataFile> rawFiles) throws Exception {

    JSONObject result = new JSONObject();

    JSONObject manifest = new JSONObject();
    JSONArray raw = new JSONArray();
    for (DataFile file : rawFiles) {
      JSONObject fileJson = new JSONObject();
      fileJson.put("filename", file.getFilename());
      fileJson.put("startDate", DateTimeUtils.toIsoDate(file.getStartDate()));
      fileJson.put("endDate", DateTimeUtils.toIsoDate(file.getEndDate()));
      raw.put(fileJson);
    }
    manifest.put("raw", raw);

    JSONArray datasetArray = new JSONArray();
    for (ExportOption exportOption : exportOptions) {
      JSONObject datasetObject = new JSONObject();
      datasetObject.put("destination", exportOption.getName());
      datasetObject.put("filename",
          dataset.getName() + exportOption.getFileExtension());
      datasetArray.put(datasetObject);
    }

    manifest.put("dataset", datasetArray);

    JSONObject metadata = DataSetDB.getMetadataJson(conn, dataset);
    Properties appConfig = ResourceManager.getInstance().getConfig();

    metadata.put("quince_information",
        "Data processed using QuinCe version " + appConfig.getProperty("version"));
    manifest.put("metadata", metadata);

    result.put("manifest", manifest);
    return result;
  }

  /**
   * Build a ZIP file containing a full dataset export, including
   * the raw files used to build the dataset and a manifest containing
   * metadata and details of the files.
   *
   * The {@code exportOption} defines the export format to be used. If
   * this is {@code null}, all formats will be exported.
   *
   * @param conn A database connection
   * @param dataset The dataset to export
   * @param exportOption The export option to use
   * @return The export ZIP file
   * @throws Exception All exceptions are propagated upwards
   */
  public static byte[] buildExportZip(Connection conn, Instrument instrument,
      DataSet dataset, ExportOption exportOption) throws Exception {

    ByteArrayOutputStream zipOut = new ByteArrayOutputStream();
    ZipOutputStream zip = new ZipOutputStream(zipOut);

    String dirRoot = dataset.getName();

    List<ExportOption> exportOptions;
    if (null != exportOption) {
      exportOptions = new ArrayList<ExportOption>(1);
      exportOptions.add(exportOption);
    } else {
      exportOptions = ExportConfig.getInstance().getOptions();
    }

    for (ExportOption option : exportOptions) {
      // Add the main dataset file
      String datasetPath = dirRoot + "/dataset/" + option.getName() +
          "/" + dataset.getName() + option.getFileExtension();

      ZipEntry datasetEntry = new ZipEntry(datasetPath);
      zip.putNextEntry(datasetEntry);
      zip.write(getDatasetExport(conn, instrument, dataset, option));
      zip.closeEntry();
    }

    List<Long> rawIds = dataset.getSourceFiles(conn);
    List<DataFile> files = DataFileDB.getDataFiles(conn,
        ResourceManager.getInstance().getConfig(), rawIds);

    for (DataFile file : files) {
      String filePath = dirRoot + "/raw/" + file.getFilename();

      ZipEntry rawEntry = new ZipEntry(filePath);
      zip.putNextEntry(rawEntry);
      zip.write(file.getBytes());
      zip.closeEntry();
    }

    // Manifest
    JSONObject manifest = makeManifest(conn, dataset, exportOptions, files);
    ZipEntry manifestEntry = new ZipEntry(dirRoot + "/manifest.json");
    zip.putNextEntry(manifestEntry);
    zip.write(manifest.toString().getBytes());
    zip.closeEntry();

    // Write the response
    zip.close();
    byte[] outBytes = zipOut.toByteArray();
    zipOut.close();

    return outBytes;
  }
}
