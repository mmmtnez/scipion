/***************************************************************************
 *
 * @author: Jesus Cuenca (jcuenca@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

import ij.gui.GenericDialog;


public class ContrastEnhancePlugin extends Plugin {
	// default values - @see xmipptomo.Plugin.collectParameters, ij.plugin.ContrastEnhancer
	boolean equalize, normalize, processStack=true, useStackHistogram;
	static double saturated = 0.5;
	
	public static String COMMAND="Enhance Contrast";
	
	@Override
	public String getCommand(){
		return COMMAND;
	}
	
	// adapt radius if the original image was resized 
	@Override
	public String getOptions(){
		return "Saturated Pixels:=" + saturated + "Normalize=" + normalize + "Equalize Histogram=" + equalize +
		"Slices=" + processStack + "Use Stack Histogram=" + useStackHistogram;
	}
	
	@Override
	public void collectParameters(GenericDialog gd) {
		if(gd != null){
			saturated=gd.getNextNumber();
			normalize=gd.getNextBoolean();
			equalize=gd.getNextBoolean();
			processStack=gd.getNextBoolean();
			useStackHistogram=gd.getNextBoolean();
		}
	}

	
	public String toString(){
		return getOptions();
	}

}
