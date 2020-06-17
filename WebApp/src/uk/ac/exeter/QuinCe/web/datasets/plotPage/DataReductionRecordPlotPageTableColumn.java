package uk.ac.exeter.QuinCe.web.datasets.plotPage;

import uk.ac.exeter.QuinCe.data.Dataset.DataReduction.DataReductionRecord;
import uk.ac.exeter.QuinCe.data.Dataset.QC.Flag;
import uk.ac.exeter.QuinCe.utils.StringUtils;

public class DataReductionRecordPlotPageTableColumn
  implements PlotPageTableColumn {

  private final DataReductionRecord record;

  private final String parameterName;

  /**
   * Create a column value for a parameter from a {@link DataReductionRecord}.
   *
   * @param record
   *          The data reduction record.
   * @param parameterName
   *          The parameter.
   */
  public DataReductionRecordPlotPageTableColumn(DataReductionRecord record,
    String parameterName) {

    this.record = record;
    this.parameterName = parameterName;
  }

  @Override
  public long getId() {
    return record.getMeasurementId() + parameterName.hashCode();
  }

  @Override
  public String getValue() {
    return String.valueOf(record.getCalculationValue(parameterName));
  }

  @Override
  public boolean getUsed() {
    return false;
  }

  @Override
  public String getQcMessage() {
    return StringUtils.collectionToDelimited(record.getQCMessages(), ";");
  }

  @Override
  public boolean getFlagNeeded() {
    return record.getQCFlag().equals(Flag.NEEDED);
  }

  @Override
  public Flag getQcFlag() {
    return record.getQCFlag();
  }
}
