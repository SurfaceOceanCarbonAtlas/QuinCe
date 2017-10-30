package uk.ac.exeter.QuinCe.data.Files;

import java.time.LocalDateTime;

import uk.ac.exeter.QuinCe.data.Instrument.FileDefinition;
import uk.ac.exeter.QuinCe.data.Instrument.FileDefinitionException;
import uk.ac.exeter.QuinCe.data.Instrument.RunTypes.RunTypeCategory;

/**
 * Class representing a specific line in a data file
 * @author Steve Jones
 *
 */
public class DataFileLine {

	/**
	 * The file that the line is from
	 */
	private DataFile file;
	
	/**
	 * The line number
	 */
	private int line;
	
	/**
	 * Basic constructor
	 * @param file The data file
	 * @param line The line number
	 */
	public DataFileLine(DataFile file, int line) {
		this.file = file;
		this.line = line;
	}
	
	/**
	 * Get the name of the file that the line is in
	 * @return The file name
	 */
	public String getFilename() {
		return file.getFilename();
	}
	
	/**
	 * Get the line number
	 * @return The line number
	 */
	public int getLine() {
		return line;
	}
	
	/**
	 * Get the date of the line
	 * @return The date
	 * @throws DataFileException If the date cannot be extracted
	 */
	public LocalDateTime getDate() throws DataFileException {
		return file.getDate(line);
	}
	
	/**
	 * Determines whether or not this line should be ignored based on its Run Type.
	 * 
	 * <p>
	 *   If the file does not have Run Types, lines should never be ignored.
	 * </p>
	 * 
	 * @return {@code true} if the line should be ignored; {@code false} if it should be used
	 * @throws DataFileException If the data cannot be extracted from the file
	 * @throws FileDefinitionException If the run type is invalid
	 */
	public boolean isIgnored() throws DataFileException, FileDefinitionException {
		
		boolean ignored = false;
		
		FileDefinition fileDefinition = file.getFileDefinition();
		
		// If a file does not have run types, it can always be used
		// for a measurement. So we only check files that have them.
		if (fileDefinition.hasRunTypes()) {
			RunTypeCategory runType = file.getRunType(line);
			ignored = runType.equals(RunTypeCategory.IGNORED_CATEGORY);
		}
			
		return ignored;
	}
}
