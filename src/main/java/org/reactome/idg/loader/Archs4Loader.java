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
import org.springframework.stereotype.Component;

// NOTE:  for loading data, start MySQL with these options:  --innodb-autoinc-lock-mode=2 --innodb_write_io_threads=16 --innodb_log_buffer_size=500000000 
// Then create schema with: "CREATE SCHEMA correlation_db DEFAULT CHARACTER SET utf8mb4 ;"

// This class can probably be modified to become a more generic data loader for other types of correlation matrices. We probably won't need
// separate Loaders for ARCHS4, Harmonizome, etc... unless we want to preserve the other identifiers (UniProt IDs, NCBI Gene IDs) that the Harmonizome files may contain.

/**
 * Loads data from ARCHS4 correlation file into MySQL.
 * @author sshorser
 *
 */
@Component
public class Archs4Loader
{
	private static final Logger logger = LogManager.getLogger();
	
	private static final String PATH_TO_TMP_FILE = "/tmp/data_for_idg";

	private static final int CHUNK_SIZE = 1000000;
	
	private String filePath;
	
	@Autowired
	GeneCorrelationDAO dao;
	
	@Autowired
	ProvenanceDAO provenanceDao;
	
	@Autowired
	CommonsPool2TargetSource daoPool;
	
	public Archs4Loader()
	{
		
	}
	
	public Archs4Loader(String filePath)
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
				String[] parts = geneSymbolHeader.split(",");
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
				
				Provenance archs4Provenance = new Provenance();
				archs4Provenance.setName("ARCHS4");
				archs4Provenance.setUrl("https://amp.pharm.mssm.edu/archs4/download.html");
				archs4Provenance.setCategory("Pairwise pearson correlation of genes across expression samples.");
				archs4Provenance.setSubcategory("Human genes");
				Provenance provenanceToUse = provenanceDao.addProvenance(archs4Provenance);
//				Provenance provenanceToUse = dao.addProvenance(archs4Provenance);
				AtomicInteger recordCount = new AtomicInteger(0);
				int lineStartOffset = 1;
				
				int numWorkers = 1;
				
				Writer out = new FileWriter(tempFileName, false);
				BufferedWriter writer = new BufferedWriter(out);
				while (scanner.hasNext())
				{
					
					String line = scanner.nextLine();
					parts = line.split(",");
					String currentGeneSymbol = parts[0].substring(1, parts[0].length()-1);

					int j = 1;
//					for (int j = 1; j <= numWorkers; j++)
					{
						// define the range of values each worker will operate on. This was more useful when there were multiple worker threads reading the file.
						int lineWidth = parts.length - lineStartOffset;
						int segmentWidth = lineWidth / numWorkers;
						final int startIndex = lineStartOffset + (j * segmentWidth) - segmentWidth;
						final int endIndex = (
												startIndex + segmentWidth // - lineStartOffset
												// This other part of the expression is to add any extra items to the last worker.
												+ (j == numWorkers ? lineWidth % numWorkers : 0)
											);
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

								if (lineCount % CHUNK_SIZE == 0)
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
				}
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
}
