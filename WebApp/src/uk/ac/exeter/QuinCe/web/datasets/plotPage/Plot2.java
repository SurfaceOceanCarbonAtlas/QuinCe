package uk.ac.exeter.QuinCe.web.datasets.plotPage;

import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import uk.ac.exeter.QuinCe.data.Dataset.QC.Flag;
import uk.ac.exeter.QuinCe.data.Instrument.FileDefinition;
import uk.ac.exeter.QuinCe.utils.DateTimeUtils;

public class Plot2 {

  /**
   * Gson generator for the main plot data
   */
  private static Gson MAIN_DATA_GSON;

  /**
   * Gson generator for the flag data
   */
  private static Gson FLAGS_GSON;

  /**
   * The source data for the plot
   */
  private final PlotPage2Data data;

  /**
   * The column ID of the X axis
   */
  private ColumnHeading xAxis = null;

  /**
   * The column ID of the Y axis
   */
  private ColumnHeading yAxis = null;

  /**
   * The plot values.
   */
  private TreeSet<PlotValue> plotValues = null;

  static {
    MAIN_DATA_GSON = new GsonBuilder()
      .registerTypeAdapter(PlotValue.class, new MainPlotValueSerializer())
      .create();

    FLAGS_GSON = new GsonBuilder()
      .registerTypeAdapter(PlotValue.class, new FlagPlotValueSerializer())
      .create();
  }

  /**
   * Constructor with minimum information.
   *
   * @param data
   *          The data for the plot.
   * @param xAxis
   *          The initial X axis ID
   * @param yAxis
   *          The initial Y axis ID
   * @throws Exception
   */
  protected Plot2(PlotPage2Data data, ColumnHeading xAxis, ColumnHeading yAxis)
    throws Exception {
    this.data = data;
    this.xAxis = xAxis;
    this.yAxis = yAxis;
  }

  /**
   * Get the x Axis label.
   *
   * @return The x axis label.
   */
  public String getXaxis() {
    return null == xAxis ? "" : xAxis.getHeading();
  }

  /**
   * Get the y Axis label.
   *
   * @return The y axis label.
   */
  public String getYaxis() {
    return null == yAxis ? "" : yAxis.getHeading();
  }

  /**
   * Get the JSON data for the main plot
   *
   * @return The main plot data.
   * @throws Exception
   */
  public String getMainData() {
    return null == plotValues ? "[]" : MAIN_DATA_GSON.toJson(plotValues);
  }

  /**
   * Get the JSON data for the flags plot
   *
   * @return The flags data
   * @throws Exception
   */
  public String getFlagData() {

    String result = "[]";

    if (null != plotValues) {
      List<PlotValue> flagValues = plotValues.stream()
        .filter(x -> x.inFlagPlot()).collect(Collectors.toList());

      result = FLAGS_GSON.toJson(flagValues);
    }

    return result;
  }

  protected void makePlotValues() throws Exception {

    TreeMap<LocalDateTime, PlotPageTableColumn> xValues = data
      .getColumnValues(xAxis);

    TreeMap<LocalDateTime, PlotPageTableColumn> yValues = data
      .getColumnValues(yAxis);

    plotValues = new TreeSet<PlotValue>();

    for (LocalDateTime time : xValues.keySet()) {
      if (yValues.containsKey(time)) {

        PlotPageTableColumn x = xValues.get(time);
        PlotPageTableColumn y = yValues.get(time);

        PlotValue plotValue;

        if (xAxis.getId() == FileDefinition.TIME_COLUMN_ID) {
          plotValue = new PlotValue(DateTimeUtils.dateToLong(time), time,
            Double.parseDouble(y.getValue()),
            y.getQcFlag().equals(Flag.FLUSHING),
            y.getFlagNeeded() ? Flag.NEEDED : y.getQcFlag());
        } else {
          plotValue = new PlotValue(DateTimeUtils.dateToLong(time),
            Double.parseDouble(x.getValue()), Double.parseDouble(y.getValue()),
            y.getQcFlag().equals(Flag.FLUSHING),
            y.getFlagNeeded() ? Flag.NEEDED : y.getQcFlag());
        }

        plotValues.add(plotValue);
      }
    }
  }

  /**
   * Initialise the plot and its data. Called from the front end when the QC
   * page is loaded.
   *
   * @throws Exception
   *           If any error occurs
   */
  public void init() {
    try {
      makePlotValues();
    } catch (Exception e) {
      data.error(e.getMessage());
    }
  }

  public String getDataLabels() {
    List<String> labels = new ArrayList<String>(4);
    labels.add(xAxis.getHeading());
    labels.add("ID");
    labels.add("GHOST");
    labels.add(yAxis.getHeading());

    return new Gson().toJson(labels);

  }

  public String getFlagLabels() {
    List<String> labels = new ArrayList<String>(4);
    labels.add(xAxis.getHeading());
    labels.add("BAD");
    labels.add("QUESTIONABLE");
    labels.add("NEEDED");

    return new Gson().toJson(labels);
  }
}