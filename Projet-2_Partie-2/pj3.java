//ce programme ecrit dans le fichier project.txt


import java.awt.*;
import java.io.*;
import java.awt.event.*;


public class pj3 extends Dialog implements MouseListener
{
    String dfname = "project.txt";
    /////////////////////////////////
    String dfname2 = "project2.txt";//contient la fonction psi
    //    String dfname3 = "project3.txt";//contient la valeur de Rr
    ////////////////////////////////
    int nv, nt;
    Triangle t[];
    Vertex v[];
    int offset = 80, hmax =380 , vmax =180, downoffset = 40 + offset;
    double scale, oldscale, xmin,xmax,ymin,ymax,oldxmin,oldxmax,oldymin
          ,oldymax,zoom0x,zoom1x,zoom0y,zoom1y,R;
    //////////////// /////////////////////
    double H=2 ,h=0.4 ;
    //////////////////////////////////////
   Button loadtriau, runjob;
    boolean mesh=false;
    TextField textfield0,textfield1,textfield2,textfield3,textfield4;
    //ajout de textfield4
    Label rhs, cx,cy, temper,val1 ;//ajout de val1
    String rhsval="0.5*exp(-(1./2)*(x*x+y*y)/(R*R))", cxval="0.5",cyval="0.6", temperval="250",Rr="0.08" ;//ajout de Rr
    Checkbox checkbox0, checkbox1;
    CheckboxGroup pbcheckboxes;
	int dwhat0=1;
	 
public class Vertex { int where; double x,y; void Vertex(){where=0; x=0; y=0;}}
public class Triangle { int where; int[] v; void Triangle(){where=0; }}
    
        
public  pj3(Frame parent) 
{
        super(parent, "Project Data Input Interface", false);
        this.setSize(hmax+2*offset, vmax+4*offset);
}

public class  getfilename extends Frame
{ 
  FileDialog fd;
  public String fname;
    public getfilename() { }
    public String getthename()
    {   
        fd = new FileDialog(this, "Choose a mesh",FileDialog.LOAD);
        fd.show();
        fname = fd.getFile();
        return fname;
    }
}


public  void readtriau(String fname){
        //int i;
        try{
        FileInputStream filename = new FileInputStream(fname);
        Reader r = new BufferedReader(new InputStreamReader(filename));
   		StreamTokenizer file = new StreamTokenizer(r);

      file.parseNumbers();
         
         file.nextToken(); nv = (int)file.nval;
         file.nextToken(); nt= (int)file.nval;
         t = new Triangle[nt];
         v = new Vertex[nv];
         
         for(int i = 0; i<nv;i++) 
         {
            v[i] = new Vertex();
            file.nextToken(); v[i].x = file.nval;
            file.nextToken(); v[i].y = file.nval;
            file.nextToken(); v[i].where = (int)file.nval;
         }
         for(int i = 0; i<nt;i++) 
         {
            t[i] = new Triangle(); t[i].v = new int[3];
            for(int j=0;j<3;j++)
            {
                file.nextToken(); 
                t[i].v[j] = (int)file.nval -1;
            }
            file.nextToken(); t[i].where = (int)file.nval;
         }
         xmin=v[0].x; xmax = v[0].x; ymin=v[0].y; ymax = v[0].y; 
         for(int i = 1; i<nv;i++) 
         {
            if(v[i].x < xmin) xmin = v[i].x;
            if(v[i].x > xmax) xmax = v[i].x;
            if(v[i].y < ymin) ymin = v[i].y;
            if(v[i].y > ymax) ymax = v[i].y;
         }
         
         scale = hmax / (xmax - xmin);  
         if( vmax/ (ymax - ymin) < scale ) scale = vmax/ (ymax - ymin);
         oldscale = scale;
         oldxmin = xmin; oldxmax = xmax;
         oldymin = ymin; oldymax = ymax;
        }
        catch (IOException Ex)
        {
            System.out.println(Ex.getMessage());
        }
    }

public void rline( Graphics g, double x0, double x1, double y0, double y1)
{
                int h0 = offset+(int)(scale * (x0 - xmin));
                int h1 = offset+(int)(scale * (x1 - xmin));
                int v0 = downoffset+(int)(scale * (y0 - ymin));
                int v1 = downoffset+(int)(scale * (y1 - ymin));
                g.drawLine(h0,v0,h1,v1);
}

public void plottriau( Graphics g)
    {            
        g.setColor(Color.blue);
         for(int k = 0; k<nt;k++)
            for(int j = 0; j<3; j++)
            { 
                int i = t[k].v[j];
                int jp = j+1; if(jp==3) jp = 0;
                int ip = t[k].v[jp];
                rline(g, v[i].x, v[ip].x, v[i].y, v[ip].y);
            }
    
    }   

public void myinterface(){
        this.setLayout(new FlowLayout());
        
        loadtriau = new Button("Load Mesh");this.add(loadtriau);
        runjob = new Button("Run Job");  	this.add(runjob);
        
        rhs = new Label("f =");
        textfield0 = new TextField(rhsval,20) ;
        this.add(rhs);  					this.add(textfield0);
        
        cx = new Label("Cx=");
        textfield1 = new TextField(cxval,20) ;
        this.add(cx); 				this.add(textfield1);
 
         cy = new Label("Cy=");
        textfield2 = new TextField(cyval,20) ;
        this.add(cy); 				this.add(textfield2);

         temper = new Label("Temperature=");
        textfield3 = new TextField(temperval,20) ;
        this.add(temper); 				this.add(textfield3);

	
	val1 = new Label("Rr=");
        textfield4 = new TextField(Rr,20) ;
        this.add(val1); 				this.add(textfield4);
	

        pbcheckboxes = new CheckboxGroup();
        checkbox0 = new Checkbox("Pb direct",pbcheckboxes,dwhat0==1);
        checkbox1 = new Checkbox("Pb inverse",pbcheckboxes,dwhat0==0);
         this.add(checkbox0);				this.add(checkbox1);

        this.pack();
        
        Graphics g = this.getGraphics();
    
        loadtriau.addActionListener(new ActionListener(){
            public void actionPerformed(ActionEvent e)
            {  
	            String name ;
	            getfilename ff = new getfilename();
	            name = ff.getthename();
	            readtriau(name);
	            mesh = true;
            	repaint();
            }
        });    
	    runjob.addActionListener(new ActionListener(){
            public void actionPerformed(ActionEvent e)
            { 
	   	    savemyparam();
		    String line;
		    Runtime rt = Runtime.getRuntime();
                    try{ 
		   
			System.out.println("exec : pj3");

			Process  ax = rt.exec("./Q8");
		
			InputStream ao = ax.getInputStream() ;
			BufferedReader br = new BufferedReader(new InputStreamReader(ao));
			while (( line= br.readLine()) != null) 
			    {   System.out.println(line);  }
		    } 
                        catch (IOException Ex)
                            { 
				System.out.println("catch execption exec ");
                                System.out.println(Ex.getMessage());
				System.exit(1);             
			    } 
				System.out.println(" Exit normal  ");
            		System.exit(0);             
              }
        }); 
        this.addMouseListener(this);
        addMouseListener(this);
}

/*  marche pas; ce serait bien de reparer, (mais pas obligatoire) */
/*
public void getrunparam()
    {   
        try{
            FileInputStream filename = new FileInputStream(dfname);
            StreamTokenizer file = new StreamTokenizer(filename);
            file.parseNumbers();             
            file.nextToken(); rhsval = file.sval;
            file.nextToken(); cxval = file.sval;
            file.nextToken(); cyval = file.sval;
             file.nextToken(); temperval = file.sval;
        }
        catch (IOException Ex)
        {
            System.out.println(Ex.getMessage());
        }
    }
*/    

public void savemyparam() 
    {   String rhsval = textfield0.getText();
        String cxval = textfield1.getText();
        String cyval = textfield2.getText();
        String temperval = textfield3.getText();
	String Rr = textfield4.getText();
        try{
            FileOutputStream filename= new FileOutputStream(dfname);
	    ////:dfname="project.txt"
	    /////////////////////////////
	    FileOutputStream filename2= new FileOutputStream(dfname2);

	    //	    FileOutputStream filename3= new FileOutputStream(dfname3);

	    ///dfname2="project2.txt
	    ///dfname3="project3.txt

            PrintStream ffile = new PrintStream(filename);
            ffile.print("\"");
            ffile.print(rhsval);
            ffile.println("\"");
            ffile.println(cxval);
            ffile.println(cyval);
	    ffile.println(temperval);
	    if(checkbox0.getState()) ffile.println("0");
	    else ffile.println("1");

	    /////////
	    ffile.println(Rr);
	    /////////

	    PrintStream ffile2 = new PrintStream(filename2);
            //ffile2.print("\"");
            ffile2.print(rhsval);
            ffile2.println("");//blanc sinon l'interpreteur repond
	    //on attendait le symbole ' '  .
	    //ffile.println(cxval);
            //ffile.println(cyval);
            //ffile.println(temperval);
            if(checkbox0.getState()) ffile.println("0");
	    else ffile.println("1");

	    //PrintStream ffile3 = new PrintStream(filename3);
            
	    //	    ffile2.print(Rr);
	    //ffile2.print("\"");
            //ffile2.print(rhsval);
            //ffile2.println("");//blanc sinon l'interpreteur repond
	    //on attendait le symbole ' '  .
	    //ffile.println(cxval);
            //ffile.println(cyval);
            //ffile.println(temperval);
            //if(checkbox0.getState()) ffile.println("0");
	    //else ffile.println("1");




         }
        catch (IOException Ex)
        {  System.out.println(Ex.getMessage());  }
    }   
    
public void paint( Graphics g ) 
{  if(mesh) plottriau(g);   }
            
public void mousePressed(MouseEvent e)
    {
	int x=e.getX(), y=e.getY();    
        zoom0x = (x-offset)/scale+xmin;
        zoom0y = (y-downoffset)/scale+ymin;
    }   

public void mouseReleased(MouseEvent e)
    {   int x=e.getX(), y=e.getY();
        zoom1x = (x-offset)/scale+xmin;
        zoom1y = (y-downoffset)/scale+ymin;
        double w0 =zoom0x-zoom1x;
        if(w0 < 0)  w0 = - w0;
        if(zoom0y-zoom1y < 0) w0 -= zoom0y-zoom1y; else w0 += zoom0y-zoom1y;
        if(w0 > 1)
        {
            xmin = zoom0x;
            xmax = zoom1x;
            ymin = zoom0y;
            ymax = zoom1y;
            scale = hmax / (xmax - xmin);
        }
        repaint();
    }       

public void mouseEntered(MouseEvent e) {;}
public void mouseExited(MouseEvent e)  {;}
public void mouseClicked(MouseEvent e) {;}

public static void main(String[] args){
        Frame f = new Frame("Java Project Interface");
        pj3 b = new pj3(f);
	// b.getrunparam();
        b.myinterface();
        b.setSize(840,300);
        b.show();
        
    }


}
