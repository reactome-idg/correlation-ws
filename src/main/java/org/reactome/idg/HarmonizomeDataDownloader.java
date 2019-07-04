package org.reactome.idg;

import java.io.IOException;
import java.net.URI;
import java.nio.file.Files;
import java.nio.file.Paths;

import org.apache.http.HttpStatus;
import org.apache.http.client.ClientProtocolException;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClientBuilder;
import org.apache.http.util.EntityUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class HarmonizomeDataDownloader
{

	private URI urlToFile;
	private String datasetName;
	private String dataCategory;
	
	private static final Logger logger = LogManager.getLogger();
	
	public HarmonizomeDataDownloader(URI url, String name, String category)
	{
		this.urlToFile = url;
		this.datasetName = name;
		this.dataCategory = category;
	}
	
	public void downloadFile() throws ClientProtocolException, IOException
	{
		try(CloseableHttpClient client = HttpClientBuilder.create().build();)
		{
			HttpGet get = new HttpGet(this.urlToFile);
			try(CloseableHttpResponse response = client.execute(get);)
			{
				if (HttpStatus.SC_OK == response.getStatusLine().getStatusCode())
				{
					byte[] b = EntityUtils.toByteArray(response.getEntity());
					String[] parts = this.urlToFile.getPath().split("/");
					String fileName = parts[parts.length-1];
					String outputPath = "/tmp/"+this.datasetName+fileName;
					Files.write(Paths.get(outputPath), b);
					logger.info("Data from {} has been downloaded to {}", this.urlToFile.toString(), outputPath);
				}
				else
				{
					logger.error("Non-200 response code ({}) was returned with message: {}", response.getStatusLine().getStatusCode(), response.getStatusLine().getReasonPhrase());
				}
			}
		}
	}

	public URI getUrlToFile()
	{
		return urlToFile;
	}

	public String getDatasetName()
	{
		return datasetName;
	}

	public String getDataCategory()
	{
		return dataCategory;
	}
}
