package junit.uk.ac.exeter.QuinCe.data.Dataset.DataReduction;

import java.sql.Connection;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import uk.ac.exeter.QuinCe.data.Dataset.DatasetSensorValues;
import uk.ac.exeter.QuinCe.data.Dataset.Measurement;
import uk.ac.exeter.QuinCe.data.Dataset.MeasurementValue;
import uk.ac.exeter.QuinCe.data.Dataset.DataReduction.DataReducer;
import uk.ac.exeter.QuinCe.data.Dataset.DataReduction.MeasurementValues;
import uk.ac.exeter.QuinCe.data.Dataset.DataReduction.ValueCalculators;
import uk.ac.exeter.QuinCe.data.Dataset.QC.Routines.RoutineException;
import uk.ac.exeter.QuinCe.data.Instrument.Instrument;
import uk.ac.exeter.QuinCe.data.Instrument.SensorDefinition.SensorType;
import uk.ac.exeter.QuinCe.data.Instrument.SensorDefinition.SensorTypeNotFoundException;
import uk.ac.exeter.QuinCe.utils.MissingParamException;

/**
 * A dummy {@link MeasurementValues} class that returns fixed values
 */
@SuppressWarnings("serial")
public class FixedMeasurementValues extends MeasurementValues {

  private Map<String, Double> sensorValues;

  /**
   * Initialise the object for a given measurement and the specified sensor
   * values.
   *
   * @param instrument
   *          The instrument that the values belong to.
   * @param measurement
   *          The measurement that the values are linked to.
   * @param searchId
   *          The search ID.
   * @param sensorValues
   *          The sensor values to use.
   */
  protected FixedMeasurementValues(Instrument instrument,
    Measurement measurement, String searchId,
    Map<String, Double> sensorValues) {

    super(instrument, measurement);
    this.sensorValues = sensorValues;
  }

  @Override
  public Double getValue(String sensorType,
    Map<String, ArrayList<Measurement>> allMeasurements,
    DatasetSensorValues allSensorValues, DataReducer reducer,
    ValueCalculators valueCalculators, Connection conn) throws Exception {

    return sensorValues.get(sensorType);
  }

  @Override
  public Double getValue(SensorType sensorType,
    Map<String, ArrayList<Measurement>> allMeasurements,
    DatasetSensorValues allSensorValues, DataReducer reducer,
    ValueCalculators valueCalculators, Connection conn) throws Exception {

    return getValue(sensorType.getName(), allMeasurements, allSensorValues,
      reducer, valueCalculators, conn);
  }

  @Override
  public void loadSensorValues(DatasetSensorValues allSensorValues,
    SensorType sensorType, boolean goodFlagsOnly)
    throws RoutineException, SensorTypeNotFoundException, MissingParamException,
    CloneNotSupportedException {

    // Do nothing
  }

  @Override
  public List<MeasurementValue> getAllValues() {
    return null;
  }

  @Override
  public void put(SensorType sensorType, MeasurementValue value) {
    // Ignore
  }
}