package org.reactome.idg.loader;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.Writer;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.time.Duration;
import java.time.LocalDateTime;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Scanner;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.reactome.idg.DataRepository;
import org.reactome.idg.dao.GeneCorrelationDAO;
import org.reactome.idg.dao.ProvenanceDAO;
import org.reactome.idg.model.Provenance;
import org.springframework.aop.target.CommonsPool2TargetSource;
import org.springframework.beans.factory.annotation.Autowired;

/**
 * Loads a gene-pair correlation matrix into a MySQL database.
 * @author sshorser
 *
 */
public class CorrelationMatrixLoader
{
	private int headerSize = 1;
	
	private String delimeter = ",";

	private static final Logger logger = LogManager.getLogger();
	
	private static final String PATH_TO_TMP_FILE = "/tmp/data_for_idg";

	@Autowired
	private int chunkSize;
	
	private String filePath;
	
	@Autowired
	private GeneCorrelationDAO dao;
	
	@Autowired
	private ProvenanceDAO provenanceDao;
	
	@Autowired
	CommonsPool2TargetSource daoPool;

	private Provenance provenance;
	
	public CorrelationMatrixLoader()
	{
		
	}
	
	public CorrelationMatrixLoader(String filePath)
	{
		this.filePath = filePath;
	}
	
	public void loadData() throws Exception
	{
		try(FileInputStream fis = new FileInputStream(this.filePath);
			Scanner scanner = new Scanner(fis);)
		{
			if (scanner.hasNext())
			{
				String geneSymbolHeader = scanner.nextLine();
				String[] parts = geneSymbolHeader.split(delimeter);
				String[] allGeneSymbols = new String[parts.length-1];
				for (int i = 1; i < parts.length; i++)
				{
					allGeneSymbols[i - 1] = parts[i].replaceAll("\"", "");
				}
				
				logger.info("{} gene symbols in the header.", allGeneSymbols.length);
				int maxPairs = (allGeneSymbols.length * (allGeneSymbols.length + 1))/2;
				logger.info("{} possible gene-pair correlations.",  maxPairs);
				
				final LocalDateTime startTime = LocalDateTime.now();
				AtomicInteger linesInFile = new AtomicInteger(0);
				AtomicInteger fileNum = new AtomicInteger(0);
				String tempFileName = PATH_TO_TMP_FILE;
				
				Provenance provenanceToUse = provenanceDao.addProvenance(this.provenance);
//				Provenance provenanceToUse = dao.addProvenance(archs4Provenance);
				AtomicInteger recordCount = new AtomicInteger(0);
				int lineStartOffset = 1;
				
				int numWorkers = 1;
				
				Writer out = new FileWriter(tempFileName, false);
				BufferedWriter writer = new BufferedWriter(out);
				
				// Read the rest of the header, but don't consume it. The gene symbols are always the first row, other rows may not 
				// be populated (such as containing "NA" for all Uniprot or NCBI identifiers)
				for (int i = 1; i < headerSize; i++)
				{
					scanner.nextLine();
				}
				// now read the actual lines.
				while (scanner.hasNext())
				{
					String line = scanner.nextLine();
					parts = line.split(delimeter);
					// if there are 3 header rows/columns, then the offset is 2 (because the arrays are 0-indexed), to skip the header and go straight to the data.
					int headerOffset = headerSize - 1;
					
					String currentGeneSymbol = parts[0].substring(1, parts[0].length()-1);

					int j = 1;
					// define the range of values each worker will operate on. This was more useful when there were multiple worker threads reading the file.
					int lineWidth = parts.length - lineStartOffset;
//					int segmentWidth = lineWidth / numWorkers;
					final int startIndex = headerOffset ;
					final int endIndex = lineWidth;
//						System.out.println("Worker #"+j + " startIndex : "+startIndex + " endIndex: "+endIndex + " Segment width: "+segmentWidth + " Line width: "+lineWidth);
					final List<String> subParts = Collections.synchronizedList(Arrays.asList(Arrays.copyOfRange(parts, startIndex, endIndex)));

					int i = -1;
					try
					{
						// Iterate through the values in the defined range.
						for (i = startIndex; i < endIndex; i++)
						{
							int lookupIndex = i - 1;
							String correlationValue = "";
							correlationValue = subParts.get(i-startIndex);
							String otherGeneSymbol = allGeneSymbols[lookupIndex];
							String key = DataRepository.generateKey(currentGeneSymbol, otherGeneSymbol);
							String keyParts[] = key.split("\\|");
							String g1 = keyParts[0];
							String g2 = keyParts[1];

							writer.write("\'"+g1+"\'\t\'"+g2+"\'\t\'"+correlationValue+"\'\t\'"+provenanceToUse.getId()+"\'\n");
							int count = recordCount.incrementAndGet();
							int lineCount = linesInFile.incrementAndGet();
							if (count % 1000000 == 0)
							{
								logger.info("{}M gene-pairs...", count/1000000);
								writer.flush();
							}

							if (lineCount % chunkSize == 0)
							{
								// Now it's time to load the file to the database.
								writer.flush();
								dao.loadGenePairsFromDataFile(tempFileName);
								Files.move(Paths.get(tempFileName), Paths.get(tempFileName + "_" + fileNum.getAndIncrement()));
								out = new FileWriter(tempFileName, false);
								writer = new BufferedWriter(out);
							}
						}
					}
					catch (ArrayIndexOutOfBoundsException e)
					{
						logger.info("ArrayOutOfBounds: startIndex: {} endIndex: {} i: {}", startIndex, endIndex, i);
						e.printStackTrace();
						writer.close();
						throw e;
					}
					catch (Exception e)
					{
						e.printStackTrace();
						writer.close();
						throw e;
					}
				}
				writer.flush();
				lineStartOffset++;

				writer.flush();
				// load the last file, probably is smaller than CHUNK_SIZE.
				dao.loadGenePairsFromDataFile(tempFileName);
				Files.move(Paths.get(tempFileName), Paths.get(tempFileName + "_" + fileNum.getAndIncrement()));
				writer.close();
				LocalDateTime endTime = LocalDateTime.now();
				logger.info("{} time spent loading the data.", Duration.between(startTime, endTime).toString());
			}
		}
	}

	public Provenance getProvenance()
	{
		return provenance;
	}

	public void setProvenance(Provenance provenance)
	{
		this.provenance = provenance;
	}

	public String getDelimeter()
	{
		return delimeter;
	}

	public void setDelimeter(String delimeter)
	{
		this.delimeter = delimeter;
	}

	public int getHeaderSize()
	{
		return headerSize;
	}

	public void setHeaderSize(int headerSize)
	{
		this.headerSize = headerSize;
	}

	public String getFilePath()
	{
		return filePath;
	}

	public void setFilePath(String filePath)
	{
		this.filePath = filePath;
	}
}
