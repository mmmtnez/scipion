package xmipp.ij.commons;

import ij.IJ;
import ij.WindowManager;
import ij.gui.ImageWindow;

import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import xmipp.ij.commons.XmippMenuBar.IJRequirement;

public class XmippImageWindow extends ImageWindow implements XmippIJWindow
{

	protected XmippMenuBar menu;

	public static void main(String[] args)
	{
		try
		{
			// openImageJ(Tool.VIEWER);
			//XmippStackWindow w = new XmippStackWindow(new ImagePlusLoader("/home/airen/hand.vol"));
			XmippImageWindow w = new XmippImageWindow(new ImagePlusLoader("/home/airen/coss/PPPIauxRS_afterRotation.xmp"));
			// IJ.open( "/home/airen/Coss/Xmipp/BPV_2/InputData/BPV_1386.mrc");

		}
		catch (Exception e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	private ImagePlusLoader ipl;


	public XmippImageWindow(ImagePlusLoader ipl)
	{
		this(ipl, ipl.getFileName());
	}


	public XmippImageWindow(ImagePlusLoader ipl, String title)
	{
		super(ipl.getImagePlus(), new XmippImageCanvas(ipl.getImagePlus()));
		this.ipl = ipl;
		setTitle(title);
		menu = new XmippMenuBar(this);
		setMenuBar(menu);		
		addWindowListener(new WindowAdapter()
		{
			public void windowClosing(WindowEvent winEvt)
			{

				System.exit(0);//temporarily
			}
		});
	}
	
	public void openMaskToolbar(){
		menu.runCommand("Masks Tool Bar", new IJRequirement[]{IJRequirement.IMAGEJ});
	}

	
	@Override
	public void loadData()
	{
		try
		{
				((XmippImageCanvas)getCanvas()).loadData(this);
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	@Override
	public void saveDataAs(String file) throws Exception
	{
		XmippImageConverter.writeImagePlus(imp, file);
	}

	@Override
	public void saveData() throws Exception
	{
		saveDataAs(imp.getTitle());
	}
	
	@Override
	public void windowClosing(WindowEvent e) {
		super.windowClosing(e);
		if(XmippIJUtil.getXmippImageJ() != null)
			XmippIJUtil.getXmippImageJ().close();
	}
	
	public ImagePlusLoader getImagePlusLoader()
	{
		return ipl;
	}


	@Override
	public boolean isVolume()
	{
		return false;
	}


	@Override
	public boolean isStack()
	{
		return false;
	}
	
	//overwriting ImageJ event to avoid switching menu
	public void windowActivated(WindowEvent e) {
//		if (IJ.isMacintosh())
//			this.setMenuBar(Menus.getMenuBar());
		if (IJ.debugMode) IJ.write(imp.getTitle() + ": Activated");
		if (!closed) {
			//ic.requestFocus();
			WindowManager.setCurrentWindow(this);
		}
	}

}// class XmippImageWindow
