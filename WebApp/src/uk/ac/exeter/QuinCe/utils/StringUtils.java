package uk.ac.exeter.QuinCe.utils;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.stream.Collectors;

/**
 * Miscellaneous string utilities
 *
 * @author Steve Jones
 *
 */
public final class StringUtils {

  /**
   * Private constructor to prevent instantiation
   */
  private StringUtils() {
    // Do nothing
  }

  /**
   * Converts a collection of values to a single string, with a specified
   * delimiter.
   *
   * <b>Note that this does not handle the case where the delimiter is found
   * within the values themselves.</b>
   *
   * @param collection
   *          The list to be converted
   * @param delimiter
   *          The delimiter to use
   * @return The converted list
   */
  public static String collectionToDelimited(Collection<?> collection,
    String delimiter) {

    String result = "";

    if (null != collection) {
      result = collection.stream().map(c -> c.toString())
        .collect(Collectors.joining(null == delimiter ? "" : delimiter.trim()));
    }

    return result;
  }

  /**
   * Converts a String containing values separated by semi-colon delimiters into
   * a list of String values
   *
   * <b>Note that this does not handle semi-colons within the values
   * themselves.</b>
   *
   * @param values
   *          The String to be converted
   * @param delimiter
   *          The delimiter
   * @return A list of String values
   */
  public static List<String> delimitedToList(String values, String delimiter) {
    List<String> result = null;

    if (null != values) {
      if (values.length() == 0) {
        result = new ArrayList<String>();
      } else {
        result = Arrays.asList(values.split(delimiter, 0));
      }
    }

    return result;
  }

  /**
   * Convert a delimited list of integers into a list of integers
   *
   * @param values
   *          The list
   * @param delimiter
   *          The delimiter
   * @return The list as integers
   */
  public static List<Integer> delimitedToIntegerList(String values,
    String delimiter) {

    List<Integer> result = null;

    if (values != null) {
      result = new ArrayList<Integer>();

      for (String item : delimitedToList(values, delimiter)) {
        result.add(Integer.parseInt(item));
      }
    }

    return result;
  }

  /**
   * Convert a delimited list of doubles (with {@code ;} separator) into a list
   * of doubles.
   *
   * @param values
   *          The delimited list.
   * @return The list of doubles.
   */
  public static List<Double> delimitedToDoubleList(String values) {
    return delimitedToDoubleList(values, ";");
  }

  /**
   * Convert a delimited list of double into a list of doubles
   *
   * @param values
   *          The list
   * @return The list as integers
   */
  public static List<Double> delimitedToDoubleList(String values,
    String delimiter) {

    List<Double> result = null;

    if (values != null) {
      List<String> stringList = delimitedToList(values, delimiter);
      result = new ArrayList<Double>(stringList.size());

      for (String item : stringList) {
        result.add(Double.parseDouble(item));
      }
    }

    return result;
  }

  /**
   * Convert a comma-separated list of numbers to a list of longs
   *
   * @param values
   *          The numbers
   * @return The longs
   */
  public static List<Long> delimitedToLongList(String values) {
    // TODO This is the preferred way of doing this. Make the other methods do
    // the same.

    List<Long> result = null;

    if (values != null) {
      String[] numberList = values.split(",");
      result = new ArrayList<Long>(numberList.length);

      for (String number : numberList) {
        result.add(Long.parseLong(number));
      }
    }

    return result;
  }

  /**
   * Extract the stack trace from an Exception (or other Throwable) as a String.
   *
   * @param e
   *          The error
   * @return The stack trace
   */
  public static String stackTraceToString(Throwable e) {
    StringWriter sw = new StringWriter();
    PrintWriter pw = new PrintWriter(sw);
    e.printStackTrace(pw);
    return sw.toString();
  }

  /**
   * Trims all items in a list of strings. A string that starts with a single
   * backslash has that backslash removed.
   *
   * @param source
   *          The strings to be converted
   * @return The converted strings
   */
  public static List<String> trimList(List<String> source) {
    return source.stream().map(s -> trimString(s)).collect(Collectors.toList());
  }

  public static List<String> trimListAndQuotes(List<String> source) {

    List<String> result = new ArrayList<String>(source.size());

    for (int i = 0; i < source.size(); i++) {
      String noQuotes = source.get(i).replaceAll("^\"|\"$", "");
      result.add(trimString(noQuotes));
    }

    return result;
  }

  private static String trimString(String value) {
    String trimmedValue = value.trim();
    if (trimmedValue.startsWith("\\")) {
      trimmedValue = trimmedValue.substring(1);
    }
    return trimmedValue;
  }

  /**
   * Convert a String-to-String lookup map into a String.
   * <p>
   * Each map entry is converted to a {@code key=value} pair. Each entry is
   * separated by a semi-colon.
   * </p>
   * <p>
   * <b>Note:</b> There is no handling of {@code =} or {@code ;} in the keys or
   * values.
   * </p>
   *
   * @param map
   *          The Map to be converted
   * @return The String representation of the Map
   */
  public static String mapToDelimited(Map<String, String> map) {

    StringBuilder result = new StringBuilder();

    int counter = 0;
    for (Map.Entry<String, String> entry : map.entrySet()) {
      counter++;
      result.append(entry.getKey());
      result.append('=');
      result.append(entry.getValue());

      if (counter < map.size()) {
        result.append(';');
      }
    }

    return result.toString();
  }

  /**
   * Convert a semi-colon-delimited list of {@code key=value} pairs into a Map.
   *
   * @param values
   *          The String
   * @return The Map
   * @throws StringFormatException
   *           If the String is not formatted correctly
   */
  public static Map<String, String> delimitedToMap(String values)
    throws StringFormatException {

    Map<String, String> result = new HashMap<String, String>();

    for (String entry : values.split(";", 0)) {

      String[] entrySplit = entry.split("=", 0);
      if (entrySplit.length != 2) {
        throw new StringFormatException("Invalid map format", entry);
      } else {
        result.put(entrySplit[0], entrySplit[1]);
      }
    }

    return result;
  }

  /**
   * Convert a Properties object into a JSON string
   *
   * @param properties
   *          The properties
   * @return The JSON string
   */
  public static String getPropertiesAsJson(Properties properties) {

    StringBuilder result = new StringBuilder();
    if (null == properties) {
      result.append("null");
    } else {

      result.append('{');

      int propCount = 0;
      for (String prop : properties.stringPropertyNames()) {
        propCount++;
        result.append('"');
        result.append(prop);
        result.append("\":\"");
        result.append(properties.getProperty(prop));
        result.append('"');

        if (propCount < properties.size()) {
          result.append(',');
        }
      }

      result.append('}');
    }

    return result.toString();
  }

  /**
   * Create a {@link Properties} object from a string
   *
   * @param propsString
   *          The properties String
   * @return The Properties object
   * @throws IOException
   *           If the string cannot be parsed
   */
  public static Properties propertiesFromString(String propsString)
    throws IOException {
    Properties result = null;

    if (null != propsString && propsString.length() > 0) {
      StringReader reader = new StringReader(propsString);
      Properties props = new Properties();
      props.load(reader);
      return props;
    }

    return result;
  }

  /**
   * Make a valid CSV String from the given text.
   *
   * This always performs three steps:
   * <ul>
   * <li>Surround the value in quotes</li>
   * <li>Any " are replaced with "", per the CSV spec</li>
   * <li>Newlines are replaced with semi-colons</li>
   * </ul>
   *
   * While these are not strictly necessary for all values, they are appropriate
   * for this application and the target audiences of exported CSV files.
   *
   * @param text
   *          The value
   * @return The CSV value
   */
  public static String makeCsvString(String text) {
    StringBuilder csv = new StringBuilder();

    if (null == text) {
      text = "";
    }

    csv.append('"');
    csv.append(text.trim().replace("\"", "\"\"").replaceAll("[\\r\\n]+", ";"));
    csv.append('"');

    return csv.toString();
  }

  /**
   * Generate a Double value from a String, handling thousands separators
   *
   * @param value
   *          The string value
   * @return The double value
   */
  public static Double doubleFromString(String value) {
    Double result = Double.NaN;
    if (null != value && value.trim().length() > 0) {
      result = Double.parseDouble(value.replaceAll(",", "").trim());
    }

    return result;
  }
}
