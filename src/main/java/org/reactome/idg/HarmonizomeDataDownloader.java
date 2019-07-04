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

public class HarmonizomeDataDownloader
{

	private URI urlToFile;
	private String datasetName;
	private String dataCategory;
	
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
					Files.write(Paths.get("/tmp/"+this.datasetName+fileName), b);
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
