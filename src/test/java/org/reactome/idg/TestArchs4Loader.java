package org.reactome.idg;

import static org.junit.Assert.fail;

import java.io.IOException;

import org.junit.Test;
import org.reactome.idg.config.AppConfig;
import org.reactome.idg.loader.Archs4Loader;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;

public class TestArchs4Loader
{

	@SuppressWarnings("static-method")
	@Test
	public void testLoadDataIT() throws IOException
	{
		try(AnnotationConfigApplicationContext context = new AnnotationConfigApplicationContext();)
		{
			context.register(AppConfig.class);
			context.refresh();
			Object o = context.getBean("sessionFactory");
			if (null == o)
			{
				throw new NullPointerException();
			}
			
			Archs4Loader loader = (Archs4Loader) context.getBean("archs4Loader");
			try
			{
				loader.loadData();
			}
			catch (Exception e)
			{
				e.printStackTrace();
				fail();
			}
		}
	}
}
