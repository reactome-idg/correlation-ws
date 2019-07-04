package org.reactome.idg;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import org.apache.http.client.ClientProtocolException;
import org.junit.Test;
import org.reactome.idg.config.AppConfig;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;

@SuppressWarnings("static-method")
public class TestHarmonizomeDownloader
{

	@Test
	public void testDownloadHarmonizomeFileIT() throws URISyntaxException, ClientProtocolException, IOException
	{
		URI uri = new URI("https://amp.pharm.mssm.edu/static/hdfs/harmonizome/data/kea/gene_similarity_matrix_cosine.txt.gz");
		HarmonizomeDataDownloader downloader = new HarmonizomeDataDownloader(uri, "KEA_Substrates_of_Kinases", "", "/media/sshorser/data/reactome/IDGFiles");
		downloader.downloadFile();
		String[] parts = uri.getPath().split("/");
		String fileName = parts[parts.length-1];
		Path path = Paths.get("/tmp/KEA_Substrates_of_Kinases"+fileName);
		// make sure file exists.
		assertTrue(Files.exists(path));
		// make sure file has some content
		assertTrue(Files.size(path) > 0);
	}
	
	@Test
	public void testListOfDownloaders()
	{
		try(AnnotationConfigApplicationContext context = new AnnotationConfigApplicationContext();)
		{
			context.register(AppConfig.class);
			context.refresh();
			@SuppressWarnings("unchecked")
			List<HarmonizomeDataDownloader> downloaders = (List<HarmonizomeDataDownloader>) context.getBean("harmonizomeDownloaders");
			assertNotNull(downloaders);
			assertTrue(downloaders.size() > 0);
			for (HarmonizomeDataDownloader downloader : downloaders)
			{
				System.out.println(downloader.getDatasetName() + "\t" + downloader.getDataCategory() + "\t" + downloader.getUrlToFile().toString());
			}
		}
	}
}
