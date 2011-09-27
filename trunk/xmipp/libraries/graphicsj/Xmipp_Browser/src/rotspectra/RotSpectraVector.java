/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package rotspectra;

import java.awt.Image;
import java.util.ArrayList;
import javax.swing.JFrame;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import xmipp.Filename;
import xmipp.MDLabel;
import xmipp.MetaData;

/**
 *
 * @author Juanjo Vega
 */
public class RotSpectraVector {

    JFreeChart chart;
    String filename;
    String block;
    double vector[];
    ArrayList<String> images;

    public RotSpectraVector(String block, String filename, double vector[]) {
        this.block = block;
        this.filename = filename;
        this.vector = vector;

        chart = createChart(block, vector);
        images = loadFilenames(block + Filename.SEPARATOR + filename);
    }

    static JFreeChart createChart(String label, double vector[]) {
        XYSeries series = new XYSeries(label);

        for (int i = 0; i < vector.length; i++) {
            series.add(i, vector[i]);
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);

        JFreeChart chart = ChartFactory.createXYLineChart(
                "", "", "",
                dataset, PlotOrientation.VERTICAL,
                true, true, false);
        chart.removeLegend();

        return chart;
    }

    static ArrayList<String> loadFilenames(String filename) {
        ArrayList<String> list = new ArrayList<String>();

        try {
            MetaData md = new MetaData(filename);

            long ids[] = md.findObjects();
            for (int i = 0; i < ids.length; i++) {
                list.add(md.getValueString(MDLabel.MDL_IMAGE, ids[i], true));
            }
        } catch (Exception ex) {
        }

        return list;
    }

    public Image getPreview(int w, int h) {
        return chart.createBufferedImage(w, h);
    }

    public JFrame getChart() {
        ChartPanel panel = new ChartPanel(chart);
        JFrame frame = new JFrame();
        frame.getContentPane().add(panel);
        frame.pack();

        return frame;
    }

    public int getNImages() {
        return images.size();
    }

    public ArrayList<String> getImagesFilenames() {
        return images;
    }

    public String getTooltipText() {
        return block;
    }
}
