package sedcond_test_diadem;

import java.io.File;
import java.util.List;

import org.krasnow.cng.data.ReadSWC;
import org.krasnow.cng.diadem.DiademMetric;
import org.krasnow.cng.domain.SwcTreeNode;

public class second { 
	public static void main(String[] args) throws Throwable { 
		ReadSWC read=new ReadSWC();
		File object=new File("C:/vaa3d/DiademMetric/ExampleGoldStandard.swc");
		File object1=new File("C:/vaa3d/DiademMetric/ExampleTest.swc");
		List<SwcTreeNode> list1=ReadSWC.convertSwcToBinaryTreeList(object);
		
		System.out.println(list1.size());
	    System.getProperty("java.classpath");
		DiademMetric object2=new DiademMetric(object,object1,10);
		
		double score=object2.getXYThreshold();
		System.out.println(score);
		
	} 
}