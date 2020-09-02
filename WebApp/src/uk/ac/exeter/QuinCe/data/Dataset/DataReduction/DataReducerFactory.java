package uk.ac.exeter.QuinCe.data.Dataset.DataReduction;

import java.sql.Connection;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.stream.Collectors;

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

  /**
   * Get the Data Reducer for a given variable and initialise it
   *
   * @param variable
   *          The variable
   * @return The Data Reducer
   * @throws DataReductionException
   *           If the reducer cannot be retreived
   */
  public static DataReducer getReducer(Connection conn, Instrument instrument,
    Variable variable, Properties variableProperties)
    throws DataReductionException {

    DataReducer reducer;

    try {
      switch (variable.getName()) {
      case "Underway Marine pCO₂": {
        reducer = new UnderwayMarinePco2Reducer(variable, variableProperties);
        break;
      }
      case "Underway Atmospheric pCO₂": {
        reducer = new UnderwayAtmosphericPco2Reducer(variable,
          variableProperties);
        break;
      }
      case "SailDrone Marine CO₂ NRT": {
        reducer = new SaildroneMarinePco2Reducer(variable, variableProperties);
        break;
      }
      case "SailDrone Atmospheric CO₂ NRT": {
        reducer = new SaildroneAtmosphericPco2Reducer(variable,
          variableProperties);
        break;
      }
      case "Soderman": {
        reducer = new NoReductionReducer(variable, variableProperties);
        break;
      }
      default: {
        throw new DataReductionException(
          "Cannot find reducer for variable " + variable.getName());
      }
      }
    } catch (Exception e) {
      throw new DataReductionException("Cannot initialise data reducer", e);
    }

    return reducer;
  }

  private static DataReducer getSkeletonReducer(Variable variable)
    throws DataReductionException {

    DataReducer reducer = null;

    switch (variable.getName()) {
    case "Underway Marine pCO₂": {
      reducer = new UnderwayMarinePco2Reducer(variable, null);
      break;
    }
    case "Underway Atmospheric pCO₂": {
      reducer = new UnderwayAtmosphericPco2Reducer(variable, null);
      break;
    }
    case "SailDrone Marine CO₂ NRT": {
      reducer = new SaildroneMarinePco2Reducer(variable, null);
      break;
    }
    case "SailDrone Atmospheric CO₂ NRT": {
      reducer = new SaildroneAtmosphericPco2Reducer(variable, null);
      break;
    }
    case "Soderman": {
      reducer = new NoReductionReducer(variable, null);
      break;
    }
    default: {
      throw new DataReductionException(
        "Cannot find reducer for variable " + variable.getName());
    }
    }

    return reducer;
  }

  public static List<CalculationParameter> getCalculationParameters(
    Variable variable, boolean includeCalculationColumns)
    throws DataReductionException {

    DataReducer reducer = getSkeletonReducer(variable);

    List<CalculationParameter> result;

    if (includeCalculationColumns) {
      result = reducer.getCalculationParameters();
    } else {
      result = reducer.getCalculationParameters().stream()
        .filter(x -> x.isResult()).collect(Collectors.toList());
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
