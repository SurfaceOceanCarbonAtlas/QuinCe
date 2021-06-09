package uk.ac.exeter.QuinCe.data.Dataset.DataReduction;

import java.math.BigDecimal;
import java.sql.Connection;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Properties;

import uk.ac.exeter.QuinCe.data.Dataset.DataSet;
import uk.ac.exeter.QuinCe.data.Dataset.Measurement;
import uk.ac.exeter.QuinCe.data.Instrument.Instrument;
import uk.ac.exeter.QuinCe.data.Instrument.Calibration.CalculationCoefficient;
import uk.ac.exeter.QuinCe.data.Instrument.Calibration.CalculationCoefficientDB;
import uk.ac.exeter.QuinCe.data.Instrument.Calibration.CalibrationSet;
import uk.ac.exeter.QuinCe.data.Instrument.SensorDefinition.SensorTypeNotFoundException;
import uk.ac.exeter.QuinCe.data.Instrument.SensorDefinition.Variable;

public class JapanCustomReducer extends DataReducer {

  private static List<CalculationParameter> calculationParameters = null;

  private BigDecimal baseSlope = null;

  private BigDecimal baseIntercept = null;

  private BigDecimal slopeAdjustment = null;

  private BigDecimal interceptAdjustment = null;

  public JapanCustomReducer(Variable variable,
    Map<String, Properties> properties) throws SensorTypeNotFoundException {
    super(variable, properties);
  }

  @Override
  public void preprocess(Connection conn, Instrument instrument,
    DataSet dataset, List<Measurement> allMeasurements) throws Exception {

    // Get the calibration slope information

    CalibrationSet coefficients = CalculationCoefficientDB.getInstance()
      .getMostRecentCalibrations(conn, instrument,
        allMeasurements.get(0).getTime());

    baseSlope = CalculationCoefficient
      .getCoefficient(coefficients, variable, "Base Slope")
      .getBigDecimalValue();

    baseIntercept = CalculationCoefficient
      .getCoefficient(coefficients, variable, "Base Intercept")
      .getBigDecimalValue();

    slopeAdjustment = CalculationCoefficient
      .getCoefficient(coefficients, variable, "Slope Adjustment")
      .getBigDecimalValue();

    interceptAdjustment = CalculationCoefficient
      .getCoefficient(coefficients, variable, "Intercept Adjustment")
      .getBigDecimalValue();
  }

  @Override
  public void doCalculation(Instrument instrument, Measurement measurement,
    DataReductionRecord record, Connection conn) throws Exception {

    Double intakeTemperature = measurement
      .getMeasurementValue("Intake Temperature").getCalculatedValue();
    Double equilibrationTemperature = measurement
      .getMeasurementValue("Equilibrator Temperature").getCalculatedValue();
    BigDecimal pCO2TEWet = new BigDecimal(measurement
      .getMeasurementValue("pCO₂ (wet at equilibration)").getCalculatedValue());
    BigDecimal index = new BigDecimal(measurement
      .getMeasurementValue("Measurement Index").getCalculatedValue());

    BigDecimal slopeAdjust = slopeAdjustment.multiply(index);
    BigDecimal interceptAdjust = interceptAdjustment.multiply(index);

    BigDecimal slope = baseSlope.add(slopeAdjust);
    BigDecimal intercept = baseIntercept.add(interceptAdjust);

    Double calibratedCO2 = pCO2TEWet.multiply(slope).add(intercept)
      .doubleValue();

    Double pCO2SST = Calculators.calcCO2AtSST(calibratedCO2,
      equilibrationTemperature, intakeTemperature);

    record.put("ΔT", Math.abs(intakeTemperature - equilibrationTemperature));
    record.put("pCO₂ SST", pCO2SST);
  }

  @Override
  protected String[] getRequiredTypeStrings() {
    return new String[] { "Intake Temperature", "Equilibrator Temperature",
      "pCO₂ (wet at equilibration)" };
  }

  @Override
  public List<CalculationParameter> getCalculationParameters() {
    if (null == calculationParameters) {
      calculationParameters = new ArrayList<CalculationParameter>(2);

      calculationParameters
        .add(new CalculationParameter(makeParameterId(0), "ΔT",
          "Water-Equilibrator Temperature Difference", "DELTAT", "°C", false));

      calculationParameters.add(new CalculationParameter(makeParameterId(1),
        "pCO₂ SST", "pCO₂ In Water", "PCO2TK02", "μatm", true));
    }

    return calculationParameters;
  }
}
