package uk.ac.exeter.QuinCe.data.Dataset.QC.Routines;

import uk.ac.exeter.QuinCe.data.Dataset.QC.Flag;

/**
 * A QC flag generated by a QC routine
 * @author Steve Jones
 *
 */
public class RoutineFlag extends Flag {

  /**
   * The routine that generated this flag
   */
  private Class<? extends Routine> routineClass;

  /**
   * The value required by the routine
   */
  private String requiredValue;

  /**
   * The actual value
   */
  private String actualValue;

  /**
   * Basic constructor
   * @param routine The routine that generated this flag
   * @param flag The flag
   */
  public RoutineFlag(Routine routine, Flag flag, String requiredValue, String actualValue) {
    super(flag);
    this.routineClass = routine.getClass();
    this.requiredValue = requiredValue;
    this.actualValue = actualValue;
  }
}
