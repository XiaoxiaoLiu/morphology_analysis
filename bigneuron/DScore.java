package DiademScore;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import javax.activation.MimetypesFileTypeMap;

import org.krasnow.cng.data.ReadSWC;
import org.krasnow.cng.diadem.DiademMetric;

public class DScore{
	public static void main(String[] args) throws Throwable { 
		boolean flag;
		int complete_mode=0;//1 indicates using gold SWC to complete with test SWC, 0 indicates using one gold SWC to complete with other gold SWC
		ArrayList<File> Gold_SwcFiles=new ArrayList<>();
		ReadSWC read=new ReadSWC();
		String main_path="C:/vaa3d/BTU_Allen/0401_gold163_all_soma_sort_strict_swc_only_0629/0401_gold163_all_soma_sort_strict_swc_only/";
		FileWriter fileWriter_Golds=new FileWriter("C:/vaa3d/BTU_Allen/0401_gold163_all_soma_sort_strict_swc_only_0629/0401_gold163_all_soma_sort_strict_swc_only/Result_Golds.txt", true);
		fileWriter_Golds.write("image_id"+"   "+"swc_file_name"+"   "+"gold_standard_swc_file_name"+"   "+"diadem_score");
		fileWriter_Golds.write("\r\n");
		FileWriter filefailed_Golds=new FileWriter("C:/vaa3d/BTU_Allen/0401_gold163_all_soma_sort_strict_swc_only_0629/0401_gold163_all_soma_sort_strict_swc_only/Failed_Golds_file.txt", true);
		filefailed_Golds.write("image_id"+"   "+"swc_file_name"+"   "+"gold_standard_swc_file_name"+"   "+"diadem_score");
		filefailed_Golds.write("\r\n");
		PrintWriter writer_ex_Golds = new PrintWriter(filefailed_Golds);
		
		FileWriter fileWriter_strict=new FileWriter("C:/vaa3d/BTU_Allen/0401_gold163_all_soma_sort_strict_swc_only_0629/0401_gold163_all_soma_sort_strict_swc_only/Result_consensus.txt", true);
		fileWriter_strict.write("image_id"+"   "+"swc_file_name"+"   "+"gold_standard_swc_file_name"+"   "+"diadem_score");
		fileWriter_strict.write("\r\n");
		FileWriter filefailed_strict=new FileWriter("C:/vaa3d/BTU_Allen/0401_gold163_all_soma_sort_strict_swc_only_0629/0401_gold163_all_soma_sort_strict_swc_only/Failed_consensus_file.txt", true);
		filefailed_strict.write("image_id"+"   "+"swc_file_name"+"   "+"gold_standard_swc_file_name"+"   "+"diadem_score");
		filefailed_strict.write("\r\n");
		PrintWriter writer_ex_strict = new PrintWriter(filefailed_strict);
		FileWriter fileWriter=new FileWriter("C:/vaa3d/BTU_Allen/0401_gold163_all_soma_sort_strict_swc_only_0629/0401_gold163_all_soma_sort_strict_swc_only/Result.txt", true);
		fileWriter.write("image_id"+"   "+"swc_file_name"+"   "+"gold_standard_swc_file_name"+"   "+"diadem_score");
		fileWriter.write("\r\n");
		FileWriter filefailed=new FileWriter("C:/vaa3d/BTU_Allen/0401_gold163_all_soma_sort_strict_swc_only_0629/0401_gold163_all_soma_sort_strict_swc_only/Failed_file.txt", true);
		filefailed.write("image_id"+"   "+"swc_file_name"+"   "+"gold_standard_swc_file_name"+"   "+"diadem_score");
		filefailed.write("\r\n");
		PrintWriter writer_ex = new PrintWriter(filefailed);
		File file = new File(main_path);
		if(file.isDirectory()){
			String[] filelist = file.list();
			System.out.println(filelist.length);
			for (int i = 0; i < filelist.length; i++) {
				Gold_SwcFiles.clear();
				File readfile = new File(main_path + filelist[i]);
				//File GoldSWC=Get_GoldSWC(main_path + filelist[i]);
				ArrayList<File> GoldSWC_Multi=Get_GoldSWC_Multi(main_path + filelist[i]);
				System.out.println(GoldSWC_Multi.size());
				for(int k=0;k<GoldSWC_Multi.size();k++){
					flag=false;
					File GoldSWC=GoldSWC_Multi.get(k);
					String GoldSWC_Name=GoldSWC.getName();
					if(!GoldSWC_Name.endsWith(".swc")){
						System.out.println("11111111111111111111111111111111111111111111111111111111111111111111111111111\n");
						continue;
					}
					Gold_SwcFiles.add(GoldSWC);
					if(GoldSWC_Name.contains("consensus")){

						flag=true;
					}
					if(complete_mode==1){
						File Test_Directory=Get_TestDir(main_path + filelist[i]);
						String[] readfilename=Test_Directory.list();
						if(Test_Directory.isDirectory()){
							for(int j=0;j<readfilename.length;j++){
								//System.out.println(Test_Directory.getAbsolutePath());
								File TestSWC=new File(Test_Directory.getAbsolutePath()+"/"+readfilename[j]);
								String SWC_name=TestSWC.getName();

								//System.out.println(TestSWC.getName());
								DiademMetric object2=new DiademMetric(TestSWC,GoldSWC,5);
								object2.setMicrons(true);
								double score=-1;
								try {
									object2.scoreReconstruction();
									score=object2.getFinalScore();
									//[image_id, swc_file_name, gold_standard_swc_file_name, diadem_score ]
									if(flag){
										Make_SpreadSheet(readfile.getName(),TestSWC.getName(),GoldSWC.getName(),score, fileWriter_strict);
									}else{
										Make_SpreadSheet(readfile.getName(),TestSWC.getName(),GoldSWC.getName(),score, fileWriter);
									}

								} catch (Exception e) {
									// TODO Auto-generated catch block
									score=1.1;
									if(flag){
										Save_failedfile(readfile.getName(),TestSWC.getName(),GoldSWC.getName(),score, filefailed_strict);
										e.printStackTrace(writer_ex_strict);
										filefailed_strict.write("\r\n");
									}else{
										Save_failedfile(readfile.getName(),TestSWC.getName(),GoldSWC.getName(),score, filefailed);
										e.printStackTrace(writer_ex);
										filefailed.write("\r\n");
									}

									e.printStackTrace();
								}

							}
						}else{
							System.out.println("something wrong!");
						}

					}else if(complete_mode==0){
						for(int ii=0;ii<Gold_SwcFiles.size();ii++){
							for(int jj=ii+1;jj<Gold_SwcFiles.size();jj++){
								File GoldSWC1=Gold_SwcFiles.get(ii);
								File GoldSWC2=Gold_SwcFiles.get(jj);
								DiademMetric object2=new DiademMetric(GoldSWC1,GoldSWC2,5);
								object2.setMicrons(true);
								double score=-1;
								try {
									object2.scoreReconstruction();
									score=object2.getFinalScore();
									//[image_id, swc_file_name, gold_standard_swc_file_name, diadem_score ]
									Make_SpreadSheet(readfile.getName(),GoldSWC1.getName(),GoldSWC2.getName(),score, fileWriter_Golds);

								} catch (Exception e) {
									// TODO Auto-generated catch block
									score=1.1;

									Save_failedfile(readfile.getName(),GoldSWC1.getName(),GoldSWC2.getName(),score, filefailed_Golds);
									e.printStackTrace(writer_ex_Golds);
									filefailed_Golds.write("\r\n");
									e.printStackTrace();
								}
							}
						}
					}

				}


			}
		}
		filefailed.flush();
		filefailed.close();
		filefailed_strict.flush();
		filefailed_strict.close();
		filefailed_Golds.flush();
		filefailed_Golds.close();
		fileWriter.flush();
		fileWriter.close();
		fileWriter_strict.flush();
		fileWriter_strict.close();
		fileWriter_Golds.flush();
		fileWriter_Golds.close();
	}
	static void Save_failedfile(String image_id, String swc_file_name, String gold_standard_swc_file_name, Double diadem_score,FileWriter fileWriter) throws Throwable{

		fileWriter.write(image_id+"   "+swc_file_name+"   "+gold_standard_swc_file_name+"   "+diadem_score);
		fileWriter.write("\r\n");
	}

	static void Make_SpreadSheet(String image_id, String swc_file_name, String gold_standard_swc_file_name, Double diadem_score,FileWriter fileWriter) throws Throwable{
		System.out.println(image_id+" "+swc_file_name+" "+gold_standard_swc_file_name+" "+diadem_score);

		fileWriter.write(image_id+"   "+swc_file_name+"   "+gold_standard_swc_file_name+"   "+diadem_score);
		fileWriter.write("\r\n");
	}

	static File Get_TestDir(String filepath){
		File result = null;
		File readfile = new File(filepath);

		if (!readfile.isDirectory()) {

		} else if (readfile.isDirectory()) {

			String[] subfilelist=readfile.list();
			for(int j=0;j<subfilelist.length;j++){
				File read_subfile=new File(filepath +"/"+subfilelist[j]);
				if(read_subfile.isDirectory()){

					result=new File(readfile.getAbsolutePath()+"/"+read_subfile.getName());
					//System.out.println("name=" + result.getAbsolutePath());
				}						
			}	
		}
		return result;
	}

	@SuppressWarnings("null")
	static ArrayList<File> Get_GoldSWC_Multi(String filepath){
		ArrayList<File> result_list = new ArrayList<>();
		File result = null;
		File readfile = new File(filepath);
		System.out.println(readfile.getAbsolutePath());

		if (!readfile.isDirectory()) {

		} else if (readfile.isDirectory()) {
			//System.out.println("absolutepath="+ readfile.getAbsolutePath());
			String[] subfilelist=readfile.list();
			for(int j=0;j<subfilelist.length;j++){
				File read_subfile=new File(filepath +"/"+subfilelist[j]);
				if(!read_subfile.isDirectory()){

					result=new File(readfile.getAbsolutePath()+"/"+read_subfile.getName());
					System.out.println(result.getAbsolutePath());
					result_list.add(result);

					//System.out.println("name=" + result.getAbsolutePath());
				}						
			}	
		}
		return result_list;

	}


	static File Get_GoldSWC(String filepath){
		File result = null;
		File readfile = new File(filepath);

		if (!readfile.isDirectory()) {

		} else if (readfile.isDirectory()) {
			//System.out.println("absolutepath="+ readfile.getAbsolutePath());
			String[] subfilelist=readfile.list();
			for(int j=0;j<subfilelist.length;j++){
				File read_subfile=new File(filepath +"/"+subfilelist[j]);
				if(!read_subfile.isDirectory()){

					result=new File(readfile.getAbsolutePath()+"/"+read_subfile.getName());
					//System.out.println("name=" + result.getAbsolutePath());
				}						
			}	
		}
		return result;
	}

}

