/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.scipion;

import xmipp.ij.commons.SocketClient;
import xmipp.ij.commons.InputFieldsMessageDialog;
import xmipp.utils.ScipionParams;
import java.awt.Color;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.SwingUtilities;
import xmipp.ij.commons.XmippApplication;
import xmipp.jni.MetaData;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippFileChooser;
import xmipp.utils.XmippResource;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.windows.GalleryJFrame;

/**
 *
 * @author airen
 */
public class ScipionGalleryJFrame extends GalleryJFrame {

    private String type;
    private String self;
    private int port;
    private JButton cmdbutton;
    private String sqlitefile;
    private JButton classcmdbutton;
    private String inputid;
    private HashMap<String, String> msgfields;
    private final String runNameKey = "Run name:";
    private String other;
    private JButton representativesbt;
    private InputFieldsMessageDialog dlg;
    private String tmpdir;
    private JButton createvolbt;
    private String setType;
    private static final String runProtCreateSubset = "run protocol ProtUserSubSet inputObject=%s sqliteFile='%s','%s' outputClassName=%s other='%s' label='%s'";
    
   

    
    
    public ScipionGalleryJFrame(ScipionGalleryData data) {
        super(data);
        readScipionParams((ScipionParams)data.parameters);
        setScipionImageIcon();
        
    }
      
    private void setScipionImageIcon()
    {
        
            Image img = XmippResource.getIcon("scipion_logo.png").getImage();
            setIconImage(img);

    }

    protected void readScipionParams(ScipionParams parameters)
    {
        try {
            self = ((ScipionGalleryData)data).getSelf();
            setType = ((ScipionMetaData)data.getMd()).getSetType();
            type = ((ScipionGalleryData)data).getScipionType() + "s";
            port = parameters.port;
            inputid = parameters.inputid;
            String filename = data.getFileName();
            tmpdir = new File(filename).getParent() + File.separator + "tmp";
            sqlitefile = data.getTmpFile("_selection");
            msgfields = new HashMap<String, String>();
            msgfields.put(runNameKey, "create subset");
            other = parameters.other;
            initComponents();
        } catch (Exception ex) {
            Logger.getLogger(ScipionGalleryJFrame.class.getName()).log(Level.SEVERE, null, ex);
            throw new IllegalArgumentException(ex.getMessage());
        }
    }
    
    protected boolean isClass2D()
    {
        return self.equals("Class2D");
    }
    
    protected void initComponents() {
        Icon icon = XmippResource.getIcon("fa-times.png");
        JButton closebt = new JButton("Close", icon);
        closebt.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent ae) {
                close();
            }
        });
        
        buttonspn.add(closebt);
        if(!XmippApplication.isScipion())
            return;
            
        if (type != null) {
            if(!data.isCTFMd())
            {
                cmdbutton = XmippWindowUtil.getScipionIconButton("Create " + type);
                cmdbutton.addActionListener(new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent ae) {
                        createSimpleSubset();
                    }
                });
                buttonspn.add(cmdbutton);
            }
            if(data.hasClasses())
            {
                classcmdbutton = XmippWindowUtil.getScipionIconButton("Create Classes");
                classcmdbutton.addActionListener(new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent ae) {
                        createSubsetFromClasses();

                    }
                });
                
                String repText = isClass2D() ? "Create Averages": "Create Volumes";
                representativesbt = XmippWindowUtil.getScipionIconButton(repText);
                representativesbt.addActionListener(new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent ae) {
                        createRepresentativesSubset();

                    }
                });
                
                buttonspn.add(representativesbt);
                buttonspn.add(classcmdbutton);
            }
            
            
            if(data.isCTFMd())
            {
                icon = XmippResource.getIcon("fa-cogs.png");
                JButton recalculatectfbt = XmippWindowUtil.getScipionIconButton("Recalculate CTFs");
                recalculatectfbt.addActionListener(new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent ae) {
                        createCTFsSubset();
                    }
                });
                recalculatectfbt.setIcon(icon);
                
                JButton ctfsubsetbt = XmippWindowUtil.getScipionIconButton("Create Micrographs");
                ctfsubsetbt.addActionListener(new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent ae) {
                        int size = ((ScipionGalleryData)data).getEnabledCount();
                        if (confirmCreate("Micrographs", size)) 
                        {
                            String command = String.format(runProtCreateSubset, 
                            inputid, sqlitefile, "", "SetOfMicrographs", other, getRunLabel());
                            createSubset(command, "Creating set ...");
                        }
                    }
                });
                buttonspn.add(ctfsubsetbt);
                buttonspn.add(recalculatectfbt);
            }
            if(self.equals("Volume") || self.equals("Class3D"))
            {
                createvolbt = XmippWindowUtil.getScipionIconButton("Create Volume");
                createvolbt.addActionListener(new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent ae) {
                        createVolume();
                    }
                });
                buttonspn.add(createvolbt);
                createvolbt.setVisible(!data.isTableMode());
            }
            pack();
            enableActions();
            jcbBlocks.addActionListener(new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent ae) {
                    enableActions();
                }
            });
        }
        
    }
    
    
    
    protected void createVolume()
    {
        msgfields.put(runNameKey, "ProtRegisterVolume");
        boolean register = confirmCreate("Are you sure you want to register volume in scipion?");
        if(register)
        {
            String command = String.format(runProtCreateSubset, 
                inputid, sqlitefile, "", setType, data.getSelVolId().toString() + ",Volume", getRunLabel());
            createSubset(command, "Creating set ...");
        }
    }
    
    public String getRunLabel()
    {
        return dlg.getFieldValue(runNameKey);
    }
    
    public String getDBPreffix()
    {
        return ((ScipionGalleryData)data).getPreffix();
    }
    
    protected void createSimpleSubset()
    {
        int size = 0;
                    
        if(data.hasClasses())
        {
            for(ScipionMetaData.EMObject emo: ((ScipionGalleryData)data).getEMObjects())
                if(emo.isEnabled() && emo.childmd != null)
                    size += emo.childmd.getEnabledCount();
        }
        else
            size = ((ScipionGalleryData)data).getEnabledCount();
        if (confirmCreate(type, size)) 
        {
            String command = String.format(runProtCreateSubset, 
                    inputid, sqlitefile, ((ScipionGalleryData)data).getPreffix(), String.format("SetOf%s", type), other, getRunLabel());
            createSubset(command, "Creating set ...");
        }
    }
    
    
    protected void createSubsetFromClasses()
    {
        int size = ((ScipionGalleryData)data).getEnabledCount();
                       
        if (confirmCreate("Classes", size)) {
            String command = String.format(runProtCreateSubset, 
                inputid, sqlitefile, "", setType , other, getRunLabel());

            createSubset(command, "Creating set ...");

        }
    }
    
    protected void createRepresentativesSubset()
    {
        int size = ((ScipionGalleryData)data).getEnabledCount();
                        
        if (confirmCreate("Representatives", size)) {
            String output = isClass2D()? "SetOfAverages,Representatives":"SetOfVolumes,Representatives";
            String command = String.format(runProtCreateSubset, 
            inputid, sqlitefile, "", output , other, getRunLabel());
            createSubset(command, "Creating set ...");

        }
    }
    
    protected void createCTFsSubset()
    {
        try {
            
            if(!data.hasRecalculateCTF())
            {
                XmippDialog.showError(ScipionGalleryJFrame.this, "There are no ctfs to recalculate");
                return;
            }

            ((ScipionGalleryData)data).overwrite(sqlitefile);
            final String[] command = new String[]{"ctf script???", inputid, sqlitefile};
            new Thread(new Runnable() {

                @Override
                public void run() {

                    try {

                        String output = XmippWindowUtil.executeCommand(command, false);
                        
                    } catch (Exception ex) {
                        throw new IllegalArgumentException(ex.getMessage());
                    }

                }
            }).start();
        close(false);                          
        } catch (Exception ex) {
            Logger.getLogger(ScipionGalleryJFrame.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public boolean confirmCreate(String output, int size)
    {
        String msg = String.format("<html>Are you sure you want to create a new set of %s with <font color=red>%s</font> %s?", output, size, (size > 1)?"elements":"element");
        if( ((ScipionGalleryData)data).getEnabledCount() == data.size())
            msg += "<br><font color=red>Note:</font> There are no disabled items to dismiss";
        return confirmCreate(msg);
    }
    
    public boolean confirmCreate(String msg)
    {
       
        dlg = new InputFieldsMessageDialog(ScipionGalleryJFrame.this, "Question", msg, msgfields);
        int create = dlg.action;
        return (create == InputFieldsMessageDialog.OK_OPTION);
    }

    public void reloadTableData(boolean changed)
    {
        super.reloadTableData(changed);
        enableActions();
    }
    
    

    protected void enableActions() {
        boolean isenabled = !data.isVolumeMode();
        Color color = isenabled ? XmippWindowUtil.firebrick : XmippWindowUtil.lightgrey;
        Color forecolor = isenabled ? Color.WHITE : Color.GRAY;
        if(cmdbutton != null)
        {
            
            cmdbutton.setVisible(isenabled);
            cmdbutton.setBackground(color);
            cmdbutton.setForeground(forecolor);
        }
        if(classcmdbutton != null)
        {
            isenabled = data.hasClasses() && !data.isVolumeMode();
            color = isenabled? XmippWindowUtil.firebrick: XmippWindowUtil.lightgrey; 
            forecolor = isenabled? Color.WHITE: Color.GRAY;
            classcmdbutton.setVisible( isenabled);
            classcmdbutton.setBackground(color);
            classcmdbutton.setForeground(forecolor);
            representativesbt.setVisible( isenabled);
            representativesbt.setBackground(color);
            representativesbt.setForeground(forecolor);
        }
    }
    

    @Override
    protected void changeView()
    {
        super.changeView();
        
        if(self.equals("Volume") || self.equals("Class3D"))
        {
            cmdbutton.setVisible(data.isTableMode());
            createvolbt.setVisible(!data.isTableMode());
        }

    }
  
    public boolean proceedWithChanges()
    {
        return true;
    }
    
   protected void createSubset(final String command, String msg) 
    {
        XmippWindowUtil.blockGUI(ScipionGalleryJFrame.this, msg);
        new Thread(new Runnable() {

            @Override
            public void run() {
                try {
                    ((ScipionGalleryData)data).overwrite(sqlitefile);
                    XmippWindowUtil.runCommand(command, port);
//                    String output = XmippWindowUtil.executeCommand(command, true);
                    XmippWindowUtil.releaseGUI(ScipionGalleryJFrame.this.getRootPane());
//                    if (output != null && !output.isEmpty()) 
//                        System.out.println(output);
                    close(false);
//                        XmippDialog.showInfo(ScipionGalleryJFrame.this, output);
//                        
//                    }

                } catch (Exception ex) {
                    ex.printStackTrace();
                    throw new IllegalArgumentException(ex.getMessage());
                }

            }
        }).start();
    }
        
   

    /**
	 * Open another metadata separataly *
	 */
    @Override
    public void openMetadata(final MetaData md)
    {
        try
        {
            SwingUtilities.invokeLater(new Runnable() {

                @Override
                public void run() {
                    new ScipionGalleryJFrame(new ScipionGalleryData(ScipionGalleryJFrame.this, data.parameters, (ScipionMetaData)md));
                }
            });
            
        }
        catch(Exception e)
        {
            XmippDialog.showError(this, e.getMessage());
        }
    }
    
    
    
    protected void initGalleryMenu() {
            menu = new ScipionGalleryMenu();
                    
    }
    
    protected class ScipionGalleryMenu extends GalleryMenu//To customize showj menu for scipion
    {
        @Override
        protected void createItems() throws Exception
        {
            super.createItems();
            addItem(FILE_LOAD_SEL, "Load selection ...");
            addItem(FILE_SAVE_SEL, "Save selection as ...", "save_as.gif");
        }

        @Override
        protected void handleActionPerformed(ActionEvent evt)
        {
            super.handleActionPerformed(evt);
            String cmd = evt.getActionCommand();
            try
            {
                    if (cmd.equals(FILE_LOAD_SEL))
                    {
                        if (fc.showOpenDialog(ScipionGalleryJFrame.this) != XmippFileChooser.CANCEL_OPTION)
                            loadSelection(fc.getSelectedPath());
                    }
                    if (cmd.equals(FILE_SAVE_SEL))
                    {
                        fc.setSelectedFile(new File(sqlitefile));
                         if (fc.showOpenDialog(ScipionGalleryJFrame.this) != XmippFileChooser.CANCEL_OPTION)
                            saveSelection(fc.getSelectedPath());
                    }
                
                
            
            }
            catch (Exception e)
            {
                    showException(e);
            }
        }

        protected void loadSelection(String path) {
            
                ((ScipionGalleryData)data).loadSelection(path);
                reloadTableData();
            
        }

        protected void saveSelection(String path) {
            try {
                ((ScipionGalleryData)data).overwrite(path);
            } catch (SQLException ex) {
                Logger.getLogger(ScipionGalleryJFrame.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
        
        

}
