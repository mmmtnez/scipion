package particlepicker;


import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.Stroke;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;

import javax.swing.SwingUtilities;

import particlepicker.training.model.AutomaticParticle;
import particlepicker.training.model.TrainingParticle;

public abstract class ParticlePickerCanvas extends ImageCanvas implements MouseWheelListener
{

	final static BasicStroke dashedst = new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, new float[] { 10.0f }, 0.0f);
	final static BasicStroke continuousst = new BasicStroke();

	public ParticlePickerCanvas(ImagePlus imp)
	{
		super(imp);
		// TODO Auto-generated constructor stub
	}

	public void moveTo(TrainingParticle p)
	{
		int width = (int) getSrcRect().getWidth();
		int height = (int) getSrcRect().getHeight();
		int x0 = p.getX() - width / 2;
		if(x0 < 0)
			x0 = 0;
		if(x0 + width > imp.getWidth())
			x0 = imp.getWidth() - width;
		int y0 = p.getY() - height / 2;
		if(y0 < 0)
			y0 = 0;
		if(y0 + height > imp.getHeight())
			y0 = imp.getHeight() - height;
		Rectangle r = new Rectangle(x0, y0, width, height);
		if (!getSrcRect().contains(r))
		{
			setSourceRect(r);
			repaint();
		}
	}
	

	
	/**
	 * Adds particle or updates its position if onpick. If ondeletepick removes
	 * particle. Considers owner for selection to the first particle containing
	 * point. Sets dragged if onpick
	 */

	public void mousePressed(MouseEvent e)
	{
		if (getFrame().getTool() != Tool.PICKER)
		{
			super.mousePressed(e);
			return;
		}
		int x = super.offScreenX(e.getX());
		int y = super.offScreenY(e.getY());

		if (SwingUtilities.isRightMouseButton(e))
		{
			setupScroll(x, y);
			return;
		}
		
	}
	
	public void mouseEntered(MouseEvent e)
	{
		if (getFrame().getTool() != Tool.PICKER)
		{
			super.mouseEntered(e);
			return;
		}
		setCursor(crosshairCursor);
	}

	public void mouseMoved(MouseEvent e)
	{
		if (getFrame().getTool() != Tool.PICKER)
		{
			super.mouseMoved(e);
			return;
		}
		setCursor(crosshairCursor);
	}
	
	
	public void mouseDragged(MouseEvent e)
	{

		if (getFrame().getTool() != Tool.PICKER)
		{
			super.mouseDragged(e);
			return;
		}
		if (SwingUtilities.isRightMouseButton(e))
		{
			scroll(e.getX(), e.getY());
			return;
		}
	}
	
	public void mouseReleased(MouseEvent e)
	{
		if (getFrame().getTool() != Tool.PICKER)
		{
			super.mouseReleased(e);
			return;
		}
		
	}
	
	public abstract void setActive(TrainingParticle p);
	
	public abstract ParticlePickerJFrame getFrame();
	
	
	protected void drawShape(Graphics2D g2, TrainingParticle p, boolean all)
	{
		Stroke previous = g2.getStroke();
		if(p instanceof AutomaticParticle)
			g2.setStroke(dashedst);
		int x0 = (int) getSrcRect().getX();
		int y0 = (int) getSrcRect().getY();
		int size = (int) (p.getFamily().getSize() * magnification);
		int radius = (int) (p.getFamily().getSize() / 2 * magnification);
		int x = (int) ((p.getX() - x0) * magnification);
		int y = (int) ((p.getY() - y0) * magnification);
		int distance = (int) (10 * magnification);

		if (getFrame().isShapeSelected(Shape.Rectangle) || all)
			g2.drawRect(x - radius, y - radius, size, size);
		if (getFrame().isShapeSelected(Shape.Circle) || all)
			g2.drawOval(x - radius, y - radius, size, size);
		if (getFrame().isShapeSelected(Shape.Center) || all)
		{
			g2.drawLine(x, y - distance, x, y + distance);
			g2.drawLine(x + distance, y, x - distance, y);
		}
		g2.setStroke(previous);
	}
	
	@Override
	public void mouseWheelMoved(MouseWheelEvent e)
	{
		int x = e.getX();
		int y = e.getY();

		int rotation = e.getWheelRotation();
		if (rotation < 0)
			zoomIn(x, y);
		else
			zoomOut(x, y);
		if (getMagnification() <= 1.0)
			imp.repaintWindow();

	}
	
	protected void drawLine(double alpha, Graphics2D g2)
	{
		
		int c1, c2, x1 = 0, x2 = 0, y1 = 0, y2 = 0;
		int pos = (int)(alpha/(Math.PI/4))%4;
		c2 = (pos == 0 || pos == 3)? imp.getHeight()/2: imp.getWidth()/2;
		c1 = (int) (c2 * Math.tan( alpha - pos * Math.PI/4));
		if(pos == 0)
		{
			x1 = imp.getWidth()/2 - c1;
			y1 = 0;
			x2 = imp.getWidth()/2 + c1;
			y2 = imp.getHeight() - 1;
		}
		else if(pos == 1)
		{
			x1 = 0;
			y1 = imp.getHeight()/2 - c1;
			x2 = imp.getWidth() - 1;
			y2 = imp.getHeight()/2 + c1;
		}
		else if(pos == 2)
		{
			x1 = 0;
			y1 = imp.getHeight()/2 + c1;
			x2 = imp.getWidth() - 1;
			y2 = imp.getHeight()/2 - c1;
		}
		else if(pos == 3)
		{
			x1 = imp.getWidth()/2 - c1;
			y1 = imp.getHeight() - 1;
			x2 = imp.getWidth()/2 + c1;
			y2 = 0;
		}
		x1 = (int)(x1 * magnification);
		x2 = (int)(x2 * magnification);
		y1 = (int)(y1 * magnification);
		y2 = (int)(y2 * magnification);
		System.out.printf("pos: %s alpha: %s x1: %s y1:%s x2:%s y2:%s\n", pos, Math.toDegrees(alpha), x1, y1, x2, y2);
		Color ccolor = g2.getColor();
		g2.setColor(Color.yellow);
		g2.drawLine(x1, y1, x2, y2);
		g2.setColor(ccolor);
	}
	
	
	public void updateMicrographData()
	{
		Micrograph m = getMicrograph();
		imp = m.getImagePlus();
		ImageWindow iw = (ImageWindow) getParent();
		iw.setImage(imp);
		iw.updateImage(imp);
		iw.setTitle(m.getName());
		
		if(!getFrame().getParticlePicker().getFilters().isEmpty())
			IJ.runMacro(getFrame().getParticlePicker().getFiltersMacro());
	}
	
	public abstract Micrograph getMicrograph();
	
	
	


}
