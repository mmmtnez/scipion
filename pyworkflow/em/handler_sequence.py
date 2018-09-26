# -*- coding: utf-8 -*-
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


# sequence related stuff
SEQ_TYPE=['aminoacids', 'nucleotides']
SEQ_TYPE_AMINOACIDS = 0
SEQ_TYPE_NUCLEOTIDES = 1
IUPAC_PROTEIN_ALPHABET = ['Extended Protein', 'Protein']
EXTENDED_PROTEIN_ALPHABET = 0
PROTEIN_ALPHABET = 1
IUPAC_NUCLEOTIDE_ALPHABET = ['Ambiguous DNA', 'Unambiguous DNA',
                             'Extended DNA', 'Ambiguous RNA',
                             'Unambiguous RNA']
EXTENDED_DNA_ALPHABET = 0
AMBIGOUS_DNA_ALPHABET = 1
UNAMBIGOUS_DNA_ALPHABET = 2
AMBIGOUS_RNA_ALPHABET = 3
UNAMBIGOUS_RNA_ALPHABET = 4


from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import Entrez, SeqIO
import urllib, sys
from Bio.SeqRecord import SeqRecord
import os


class SequenceHandler():
    def __init__(self, sequence=None,
                 iUPACAlphabet=0,
                 isAminoacid=True):

        self.isAminoacid = isAminoacid
        if isAminoacid:
            if iUPACAlphabet == EXTENDED_PROTEIN_ALPHABET:
                self.alphabet = IUPAC.ExtendedIUPACProtein
            else:
                 self.alphabet = IUPAC.IUPACProtein
        else:
            if iUPACAlphabet == EXTENDED_DNA_ALPHABET:
                self.alphabet = IUPAC.ExtendedIUPACDNA
            elif iUPACAlphabet == AMBIGOUS_DNA_ALPHABET:
                self.alphabet = IUPAC.IUPACAmbiguousDNA
            elif iUPACAlphabet == UNAMBIGOUS_DNA_ALPHABET:
                self.alphabet = IUPAC.IUPACUnambiguousDNA
            elif iUPACAlphabet == AMBIGOUS_RNA_ALPHABET:
                self.alphabet = IUPAC.IUPACAmbiguousRNA
            else:
                self.alphabet = IUPAC.IUPACUnambiguousRNA
        if sequence is not None:
            sequence = self._cleanSequence(sequence)
            self._sequence = Seq(sequence, alphabet=self.alphabet)
        else:
            self._sequence = None
            # type(self._sequence):  <class 'Bio.Seq.Seq'>

    def _cleanSequence(self, sequence):
        alphabet = self.alphabet.letters
        str_list = []
        for item in sequence.upper():
            if item in alphabet:
                str_list.append(item)
        return ''.join(str_list)

    def saveFile(self, fileName, seqID, sequence=None, name=None,
                 seqDescription=None, type="fasta"):
        if sequence is not None:
            self._sequence = sequence
        record = SeqRecord(self._sequence, id=seqID, name=name,
                           description=seqDescription)
        # type(record): < class 'Bio.SeqRecord.SeqRecord'>
        with open(fileName, "w") as output_handle:
            SeqIO.write(record, output_handle, type)

    def downloadSeqFromFile(self, fileName, type="fasta"):
        record = next(SeqIO.parse(fileName, type))
        return record

    def downloadSeqFromDatabase(self, seqID):
        # see http://biopython.org/DIST/docs/api/Bio.SeqIO-module.html
        # for format/databses
        print "Conneting to dabase..."
        sys.stdout.flush()
        counter=1
        retries = 5
        while counter <= retries:  # retry up to 5 times if server busy
            try:
                if self.isAminoacid:
                    url = "http://www.uniprot.org/uniprot/%s.xml"
                    format = "uniprot-xml"
                    handle = urllib.urlopen(url % seqID)
                else:
                    Entrez.email = os.environ.get('SCIPION_MAIL', None)
                    handle = Entrez.efetch(db="nucleotide", id=seqID,
                                           rettype="gb", retmode="text")
                    format = "gb"
                record = SeqIO.read(handle, format)
                break
            except Exception, e:
                counter += 1
                print "Can not connect to database retry %d of %d" % (counter, retries)
                print str(e)
                if counter == retries:
                    self.protocol.setAborted()
        return record
