package org.reactome.idg.loader;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.time.Duration;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Scanner;
import java.util.concurrent.Callable;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.atomic.AtomicInteger;

import org.reactome.idg.DataRepository;
import org.reactome.idg.dao.GeneCorrelationDAO;
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
	private static final String PATH_TO_TMP_FILE = "/tmp/data_for_idg";

	private static final int CHUNK_SIZE = 1000000;
	
	private String filePath;
	
	@Autowired
	GeneCorrelationDAO dao;
	
	@Autowired
	CommonsPool2TargetSource daoPool;
	
	public Archs4Loader()
	{
		
	}
	
	public Archs4Loader(String filePath)
	{
		this.filePath = filePath;
	}
	
//	private void changeToNextFile(int fileNum, Writer out, BufferedWriter writer) throws IOException
//	{
//		Files.move(Paths.get(PATH_TO_TMP_FILE), Paths.get(PATH_TO_TMP_FILE + "_" + fileNum));
////		writer.close();
////		out.close();
//		out = new FileWriter(PATH_TO_TMP_FILE, false);
//		writer = new BufferedWriter(out);
//	}
	
	public void loadData() throws Exception
	{
		
//		Writer out = new FileWriter(PATH_TO_TMP_FILE);
		
		
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
					allGeneSymbols[i - 1] = parts[i].replaceAll("\"", ""); //substring(1, parts[i].length() - 1);
				}

				System.out.println(allGeneSymbols.length + " gene symbols in the header.");
				int maxPairs = (allGeneSymbols.length * (allGeneSymbols.length + 1))/2;
				System.out.println( maxPairs + " possible gene-pair correlations.");
				
				final LocalDateTime startTime = LocalDateTime.now();
//				LocalDateTime cycleStartTime = startTime;
//				AtomicInteger cycleCount = new AtomicInteger(0);
				AtomicInteger linesInFile = new AtomicInteger(0);
				AtomicInteger fileNum = new AtomicInteger(0);
				String tempFileName = PATH_TO_TMP_FILE;
				
				Provenance archs4Provenance = new Provenance();
				archs4Provenance.setName("ARCHS4");
				archs4Provenance.setUrl("https://amp.pharm.mssm.edu/archs4/download.html");
				archs4Provenance.setCategory("Pairwise pearson correlation of genes across expression samples.");
				archs4Provenance.setSubcategory("Human genes");
				
//				GeneCorrelationDAO dao = new GeneCorrelationDAOImpl();
				
				Provenance provenanceToUse = dao.addProvenance(archs4Provenance);
//				dao.setCurrentProvenance(provenanceToUse);
				AtomicInteger recordCount = new AtomicInteger(0);
				int lineStartOffset = 1;
				
				int numWorkers = 1;
				
				ForkJoinPool pool = new ForkJoinPool(numWorkers);
				Collection<Callable<?>> tasks = new ArrayList<>(numWorkers);
				Writer out = new FileWriter(tempFileName, false);
				BufferedWriter writer = new BufferedWriter(out);
				while (scanner.hasNext())
				{
//					if (lineStartOffset%10 == 0)
//					{
//						System.out.println("lineOffset: "+lineStartOffset);
//					}
					
					String line = scanner.nextLine();
					parts = line.split(",");
					String currentGeneSymbol = parts[0].substring(1, parts[0].length()-1);

					
//					if (Duration.between(cycleStartTime, LocalDateTime.now()).toMinutes() >= 3)
//					{
//						System.out.println("Resting for a few moments, letting the database catch up.");
//						Thread.sleep(30 * 1000);
//						cycleStartTime = LocalDateTime.now();
//						cycleCount.set(0);
//					}
//					final LocalDateTime rateStartTime = cycleStartTime;

					
					for (int j = 1; j <= numWorkers; j++)
					{
						// define the range of values each worker will operate on.
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
//						Callable<Void> task = new Callable<Void>()
//						{
//							long localLineCount = 0 ; 
							
//							@Override
//							public Void call() throws Exception
//							{
								GeneCorrelationDAO dao1;
								int i = -1;
								try
								{
									dao1 = (GeneCorrelationDAO) daoPool.getTarget();
									dao1.setCurrentProvenance(provenanceToUse);
									dao1.setBatchSize(100000);
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
										// if we're at the end, set the batchsize to 1 to trigger commit of the remainder.
										if (i == endIndex - 1)
										{
											dao1.setBatchSize(1);
										}
//										dao1.addGenePair(g1 , g2, Double.valueOf(correlationValue));
//										synchronized(writer)
										{
//											localLineCount ++;
											writer.write("\'"+g1+"\'\t\'"+g2+"\'\t\'"+correlationValue+"\'\t\'"+provenanceToUse.getId()+"\'\n");
										}
										
//										synchronized (recordCount)
										{
											int count = recordCount.incrementAndGet();
											int lineCount = linesInFile.incrementAndGet();
//											int cCount = cycleCount.incrementAndGet();
											if (count % 1000000 == 0)
											{
												System.out.println(count/1000000 + "M gene-pairs loaded. Accumulative rate (gene-pairs / second): " + (count * 1.0) / (Duration.between(startTime, LocalDateTime.now()).toMillis()/1000) );
												writer.flush();
											}
											if (lineCount % CHUNK_SIZE == 0 || maxPairs - lineCount < CHUNK_SIZE)
											{
												// Now it's time to load the file to the database.
												writer.flush();
//												dao.loadGenePairsFromDataFile(PATH_TO_TMP_FILE);
//												changeToNextFile(fileNum.getAndIncrement(), out, writer);
												dao.loadGenePairsFromDataFile(tempFileName);
												Files.move(Paths.get(tempFileName), Paths.get(tempFileName + "_" + fileNum.getAndIncrement()));
												out = new FileWriter(tempFileName, false);
												writer = new BufferedWriter(out);
											}
										}
									}
//									System.out.println("I wrote "+localLineCount+" lines");
									daoPool.releaseTarget(dao1);
//									return null;
								}
								catch (ArrayIndexOutOfBoundsException e)
								{
									System.out.println("ArrayOutOfBounds: startIndex: "+startIndex+" endIndex: "+endIndex+" i: "+i);
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
								
//							}
//						};
//						tasks.add(task);
					}
					pool.invokeAll(tasks);
					writer.flush();
					// this removes all old tasks
					tasks.clear();
					/*
					for (int i = lineStartOffset; i < parts.length; i++)
					{
						int lookupIndex = i - 1;
						String correlationValue = "";
						correlationValue = parts[i];
						
						// Get the other gene symbol for this position in the line.
						String otherGeneSymbol = allGeneSymbols[lookupIndex];
						String key = DataRepository.generateKey(currentGeneSymbol, otherGeneSymbol);
//						final Double correlationAsDouble = Double.valueOf(correlationValue);
//						if (!correlationValues.containsKey(key))
//						if (!jedis.exists(key))
//						{
							String keyParts[] = key.split("\\|");
							String g1 = keyParts[0];
							String g2 = keyParts[1];
//							try
//							{
//								dao.addGenePair(g1 , g2, Double.valueOf(correlationValue));
								
								GeneCorrelationDAO dao1 = (GeneCorrelationDAO) daoPool.getTarget();
								dao1.setCurrentProvenance(provenanceToUse);
								dao1.addGenePair(g1 , g2, Double.valueOf(correlationValue));
//							}
//							catch (ConstraintViolationException e)
//							{
//								// no-op
//							}
							
//							correlationValues.put(key, correlationAsDouble);
//							jedis.set(key, correlationValue);
							c++;
							if (c % 10000 == 0 )
							{
//								System.out.println(jedis.dbSize() + " gene-pairs loaded.");
								System.out.println(c + " gene-pairs loaded.");
							}
//						}
						// If by some chance this key has already been seen (possibly because of a duplicated identifier), append the correlation value.
//						else 
//						{
//							System.out.println("Whoa! The GeneSymbol key " + key + " has already been seen, with value: " + jedis.get(key) );
//						}
					}
					*/
					lineStartOffset++;
				}
				writer.close();

				LocalDateTime endTime = LocalDateTime.now();
				// Now that we're done, print out the elapsed time and number of keys loaded.
//				System.out.println(correlationValues.keySet().size() + " gene-pair keys (Gene Symbol) have been loaded.");
				System.out.println(Duration.between(startTime, endTime).toString() + " time spent loading the data.");
				
				// Now, load the file into MySQL
//				dao.loadGenePairsFromDataFile(PATH_TO_TMP_FILE);
			}
		}
	}	
}
