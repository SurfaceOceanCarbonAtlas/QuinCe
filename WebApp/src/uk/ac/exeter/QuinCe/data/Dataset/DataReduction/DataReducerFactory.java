package uk.ac.exeter.QuinCe.data.Dataset.DataReduction;

import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;

import uk.ac.exeter.QuinCe.data.Instrument.Instrument;
import uk.ac.exeter.QuinCe.data.Instrument.InstrumentException;
import uk.ac.exeter.QuinCe.data.Instrument.SensorDefinition.Variable;

/**
 * Factory class for Data Reducers
 *
 * @author Steve Jones
 *
 */
public class DataReducerFactory {

  private static Map<String, Class<? extends DataReducer>> reducers;

  static {
    reducers = new HashMap<String, Class<? extends DataReducer>>();
    reducers.put("CONTROS pCO₂", ControsPco2Reducer.class);
    reducers.put("SailDrone Atmospheric CO₂ NRT",
      SaildroneAtmosphericPco2Reducer.class);
    reducers.put("SailDrone Marine CO₂ NRT", SaildroneMarinePco2Reducer.class);
    reducers.put("Soderman", NoReductionReducer.class);
    reducers.put("Underway Atmospheric pCO₂",
      UnderwayAtmosphericPco2Reducer.class);
    reducers.put("Underway Marine pCO₂", UnderwayMarinePco2Reducer.class);
  }

  /**
   * Get the Data Reducer for a given variable and initialise it
   *
   * @param variable
   *          The variable
   * @return The Data Reducer
   * @throws DataReductionException
   *           If the reducer cannot be retreived
   */
  public static DataReducer getReducer(Variable variable,
    Map<String, Properties> properties) throws DataReductionException {

    try {
      Class<? extends DataReducer> clazz = getReducerClass(variable.getName());

      Constructor<? extends DataReducer> constructor = clazz
        .getConstructor(Variable.class, Map.class);

      return constructor.newInstance(variable, properties);
    } catch (Exception e) {
      throw new DataReductionException(
        "Cannot get reducer for variable '" + variable.getName() + "'", e);
    }
  }

  private static DataReducer getSkeletonReducer(Variable variable)
    throws DataReductionException {

    return getReducer(variable, null);
  }

  public static Class<? extends DataReducer> getReducerClass(String variable)
    throws DataReductionException {

    Class<? extends DataReducer> result = reducers.get(variable);

    if (null == result) {
      throw new DataReductionException(
        "Cannot find reducer for variable " + variable);
    }

    return result;
  }

  public static List<CalculationParameter> getCalculationParameters(
    Variable variable, boolean includeCalculationColumns)
    throws DataReductionException {

    DataReducer reducer = getSkeletonReducer(variable);

    List<CalculationParameter> allParameters = reducer
      .getCalculationParameters();

    List<CalculationParameter> result = new ArrayList<CalculationParameter>(
      allParameters.size());

    for (CalculationParameter param : allParameters) {
      boolean use = true;

      if (!includeCalculationColumns && !param.isResult()) {
        use = false;
      }

      if (use) {
        result.add(param);
      }
    }

    return result;
  }

  public static Map<Variable, List<CalculationParameter>> getCalculationParameters(
    Collection<Variable> variables) throws DataReductionException {

    Map<Variable, List<CalculationParameter>> result = new HashMap<Variable, List<CalculationParameter>>();

    for (Variable variable : variables) {
      DataReducer reducer = getSkeletonReducer(variable);
      result.put(variable, reducer.getCalculationParameters());
    }

    return result;
  }

  protected static long makeParameterId(Variable variable, int sequence) {
    return variable.getId() * 10000 + sequence;
  }

  public static Variable getVariable(Instrument instrument, long parameterId)
    throws InstrumentException {

    return instrument.getVariable(parameterId / 10000);
  }

  public static CalculationParameter getVariableParameter(Variable variable,
    long parameterId) throws DataReductionException {

    int parameterIndex = (int) (parameterId % 10000);
    return getSkeletonReducer(variable).getCalculationParameters()
      .get(parameterIndex);
  }
}
