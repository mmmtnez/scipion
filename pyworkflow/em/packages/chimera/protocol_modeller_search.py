# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os

from pyworkflow import VERSION_1_2
from pyworkflow.em import PdbFile
from pyworkflow.em import Volume
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.data import Transform
from pyworkflow.em.headers import Ccp4Header
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.em.viewers.chimera_utils import \
    createCoordinateAxisFile, getProgram, runChimeraProgram, \
    chimeraPdbTemplateFileName, chimeraMapTemplateFileName, \
    chimeraScriptFileName, sessionFile
from pyworkflow.protocol.params import MultiPointerParam, PointerParam, \
    StringParam
from pyworkflow.utils.properties import Message
from pyworkflow.em.packages.chimera.protocol_base import ChimeraProtBase


class ChimeraModelFromTemplate(ChimeraProtBase):
    """Protocol to look for structures of homologous sequences of the input
        sequence using Chimera.
        Execute command *scipionwrite [model #n]* from command line in order
        to transfer the selected
        pdb to scipion. Default value is model=#0,
        model refers to the pdb file."""
    _label = 'model from template'
    _program = ""
    _version = VERSION_1_2

    # --------------------------- DEFINE param functions --------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('pdbFileTemplate', PointerParam,
                      pointerClass="PdbFile",
                      label='PDBx/mmCIF file template',
                      help="PDBx/mmCIF file template used as basic atomic "
                           "structure to model your specific sequence.")
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      condition='False',
                      label='Input Volume', allowsNull=True,
                      help="Volume to process")
        form.addSection(label='Help')
        form.addLine('''Step 1: In Chimera main menu, select Tools -> Sequence 
        -> Sequence; Select the specific chain as template to model your 
        sequence and press Show. The sequence of the selected chain will be 
        showed in an individual window.\n\nStep 2: In the sequence window menu, 
        select Edit -> Add sequence. A new window will be opened in which you 
        can paste the specific sequence that you want to model (Plain text) 
        or download it from a text file (From File), or download it from 
        UniProt by using the UniProt name/ID (From UniProt). Once selected 
        the sequence, press OK. This last sequence will appear aligned to the 
        template's sequence.\n\nStep 3: In the sequence window menu, select 
        Structure -> Modeller (homology)...; A new window for Comparative 
        Modeling with Modeller will appear. Select your specific sequence as 
        the sequence to be modeled, and the input atomic structure used as 
        template for modeling. Select Run Modeller via web service and write 
        the Modeller license key supplied. Finally, press OK.\nWAITING TIME: 
        (you may see the status of your job in chimera main window, 
        lower left corner.)\n\nStep 4: When finished click quit, 5 models 
        will appear. Choose the one you like the best. Save it as first guess 
        in Scipion by executing the Chimera command *scipionwrite [model #n]* 
        (Chimera main menu: Favorites -> Command Line).\n 
        When you use the command line scipionwrite, the Chimera session will 
        be saved by default. Additionally, you can save the Chimera session 
        whenever you want by executing the command *scipionss". You will be 
        able to restore the saved session by using the protocol chimera restore 
        session (SCIPION menu: Tools/Calculators/chimera restore session). ''')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runChimeraStep')
        self._insertFunctionStep('createOutput')



