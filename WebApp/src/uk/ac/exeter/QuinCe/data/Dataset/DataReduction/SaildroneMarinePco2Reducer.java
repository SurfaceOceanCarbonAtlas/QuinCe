package uk.ac.exeter.QuinCe.data.Dataset.DataReduction;

import java.sql.Connection;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Properties;

import uk.ac.exeter.QuinCe.data.Dataset.Measurement;
import uk.ac.exeter.QuinCe.data.Instrument.Instrument;
import uk.ac.exeter.QuinCe.data.Instrument.SensorDefinition.Variable;

/**
 * Data Reduction class for NRT Marine fCO₂ from SailDrones.
 *
 * <p>
 * Calculations from Sutton et al. 2014 (doi: 10.5194/essd-6-353-2014).
 * </p>
 *
 * @author Steve Jones
 *
 */
public class SaildroneMarinePco2Reducer extends DataReducer {

  private static List<CalculationParameter> calculationParameters = null;

  public SaildroneMarinePco2Reducer(Variable variable,
    Map<String, Properties> properties) {

    super(variable, properties);
  }

  @Override
  public void doCalculation(Instrument instrument, Measurement measurement,
    DataReductionRecord record, Connection conn) throws Exception {

    Double intakeTemperature = measurement
      .getMeasurementValue("Intake Temperature").getCalculatedValue();
    Double salinity = measurement.getMeasurementValue("Salinity")
      .getCalculatedValue();
    Double licorPressure = measurement
      .getMeasurementValue("LICOR Pressure (Equilibrator)")
      .getCalculatedValue();
    Double co2InGas = measurement
      .getMeasurementValue("xCO₂ water (dry, no standards)")
      .getCalculatedValue();

    Double pH2O = Calculators.calcPH2O(salinity, intakeTemperature);
    Double pCO2 = Calculators.calcpCO2TEWet(co2InGas, licorPressure, pH2O);
    Double fCO2 = Calculators.calcfCO2(pCO2, co2InGas, licorPressure,
      intakeTemperature);

    record.put("pCO₂", pCO2);
    record.put("fCO₂", fCO2);
  }

  @Override
  protected String[] getRequiredTypeStrings() {
    return new String[] { "Intake Temperature", "Salinity",
      "LICOR Pressure (Equilibrator)", "xCO₂ water (dry, no standards)" };
  }

  @Override
  public List<CalculationParameter> getCalculationParameters() {
    if (null == calculationParameters) {
      calculationParameters = new ArrayList<CalculationParameter>(3);
      calculationParameters
        .add(new CalculationParameter(makeParameterId(0), "pH₂O",
          "Marine Water Vapour Pressure", "RH2OX0EQ", "hPa", false, false));
      calculationParameters.add(new CalculationParameter(makeParameterId(1),
        "pCO₂", "pCO₂ In Water", "PCO2TK02", "μatm", true, false));
      calculationParameters.add(new CalculationParameter(makeParameterId(2),
        "fCO₂", "fCO₂ In Water", "FCO2XXXX", "μatm", true, false));
    }
    return calculationParameters;
  }
}
