import java.io.*;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.*;
import java.util.regex.*;
public class SXKeywordSearch 
{
       static HashMap hm = new HashMap();
       static HashMap snpmap = new HashMap(); 
       static Set<String> rec = new HashSet<String>();
       String GetRecord(String Filename) throws FileNotFoundException 
       {
            
            String Fname[] = Filename.split(":");
                 String dirname = "";
                 String[] s= null;
                 s = Fname[0].split("\\.");
                 dirname += "/";
                 dirname += s[0];
                 dirname += "/";
                 dirname += s[1];
                 dirname += "/";
                 dirname += s[2];
                 dirname += "/";
                 dirname += Fname[0];                
            RandomAccessFile rand = new RandomAccessFile("../../species/" + dirname,"rw");
            String record = new String();
            try 
            {
			rand.seek(Long.parseLong(Fname[1]));
                        record = rand.readLine();
			rand.close();
            }
	    catch (IOException e) 
            {
			e.printStackTrace();
	    }	 
            return record;
      }
     void SearchKeyword(String keyword)
     { 
          if(snpmap.containsKey( keyword.toLowerCase() ))
          {
                   try
                   {
                        rec.add(GetRecord((String)snpmap.get(keyword.toLowerCase())));
                   }
                   catch (FileNotFoundException e)
                   {
                        e.printStackTrace();
                   }
          }
          else
          {
            if(hm.containsKey(keyword.toLowerCase() ) ) 
            {
               List list = (List) hm.get(keyword.toLowerCase());
               Iterator<String> iterator = list.iterator();
	       while ( iterator.hasNext() )
               {
                   try
                   {
                        rec.add(GetRecord(iterator.next()));
                   } 
                   catch (FileNotFoundException e) 
                   {
                        e.printStackTrace();
                   }
               }
            }
            else
            {
                   Set set = hm.entrySet();
                   Iterator i = set.iterator();
                   while(i.hasNext())
                   {
                      Map.Entry me = (Map.Entry)i.next();
                      String h = (String)me.getKey();
                      List list = (List) hm.get(h);
                      Iterator<String> iterator = list.iterator();
                      Pattern p = Pattern.compile(keyword.toLowerCase());
                      Matcher m = p.matcher(h);
                      
                        if(m.find())
                        {
                                 List list1 = (List) hm.get(h.toLowerCase());
                                 Iterator<String> iterator1 = list.iterator();
                                 while ( iterator1.hasNext() )
                                {
                                    try
                                    {
                                       rec.add(GetRecord(iterator1.next()));
                                    }
                                    catch (FileNotFoundException e)
                                    {
                                       e.printStackTrace();
                                    }
                               
                         }
                        }
                    }
              }
         }
        }
       
    void LoadMap(String Filename ) throws FileNotFoundException 
    {
           FileReader fr = new FileReader(Filename);
	   BufferedReader br = new BufferedReader(fr);
	   String s;
	   try {
		while((s = br.readLine()) != null) 
		   { 
		      String filepath[] = s.split("\t");
		          if ( ! hm.containsKey( filepath[0] ) ) 
		          {
		             List list = new ArrayList( );
		             list.add(filepath[1]);
		             hm.put( filepath[0], list);
		          }
		          else 
		          {
		            List list = (List) hm.get(filepath[0]);
		            list.add( filepath[1]);
		          }
		   }
	} 
        catch (IOException e) 
        {
		e.printStackTrace();
	}
	   try {
		fr.close();
	} catch (IOException e) {
		e.printStackTrace();
	} 
    }

    void LoadSnpMap(String Filename ) throws FileNotFoundException
    {
           FileReader snpindexfile = new FileReader(Filename);
           BufferedReader br = new BufferedReader(snpindexfile);
           String s;
           try {
                while((s = br.readLine()) != null)
                   {
                      String filepath[] = s.split("\t");
                      snpmap.put( filepath[0], filepath[1]);
                   }
              } catch (IOException e) {
                e.printStackTrace();
        }
           try {
                snpindexfile.close();
        } catch (IOException e) {
                e.printStackTrace();
        }
    }

    public static void main(String args[]) throws IOException
    {
           SXKeywordSearch sobj = new SXKeywordSearch();
           sobj.LoadMap("/data/project/0811_gbrowser/data_repository/main/KeywordSearch_java/Indexfile.txt");                         
           sobj.LoadSnpMap("/data/project/0811_gbrowser/data_repository/main/KeywordSearch_java/SNPIndex.txt");  
           for(int i = 0; i < args.length; i++)
           if(!args[i].contains("+"))   
	   sobj.SearchKeyword(args[i]);
      
           Iterator iter = rec.iterator();
           while (iter.hasNext()) 
           {
                  String recvalue = (String)iter.next();
                  if(args.length == 1)
                      System.out.println(recvalue);
                  else if(args.length == 2 && (recvalue.indexOf(args[0]) != -1 || recvalue.indexOf(args[1]) != -1 ))
                      System.out.println(recvalue);
                  else if(args.length == 3 && args[1].contains("+") && recvalue.indexOf(args[0]) != -1 && recvalue.indexOf(args[1]) != -1 )
                      System.out.println(recvalue);
           } 
       
     }	 
}

