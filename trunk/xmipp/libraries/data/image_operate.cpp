/***************************************************************************
 *
 * Authors: J.M. de la Rosa Trevin   (jmdelarosa@cnb.csic.es) x 0.95
 *          Joaquin Oton             (joton@cnb.csic.es)      x 0.05
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

#include "image_operate.h"
#include "metadata_extension.h"

void minus(Image<double> &op1, const Image<double> &op2)
{
    op1() -= op2();
}

void plus(Image<double> &op1, const Image<double> &op2)
{
    op1() += op2();
}

void mult(Image<double> &op1, const Image<double> &op2)
{
    op1() *= op2();
}

void divide(Image<double> &op1, const Image<double> &op2)
{
    op1() /= op2();
}

void min(Image<double> &op1, const Image<double> &op2)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op1())
    {
        dAi(op1(), n) = XMIPP_MIN(dAi(op1(), n), dAi(op2(), n));
    }
}

void max(Image<double> &op1, const Image<double> &op2)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op1())
    {
        dAi(op1(), n) = XMIPP_MAX(dAi(op1(), n), dAi(op2(), n));
    }
}

void compare(Image<double> &op1, const Image<double> &op2)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op1())
    {
        dAi(op1(), n) = dAi(op1(), n) == dAi(op2(), n) ? 0 : (dAi(op1(), n) < dAi(op2(), n) ? -1 : 1 );
    }
}

///Be careful with integer images for relational operations...due to double comparisons
void eq(Image<double> &op1, const Image<double> &op2)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op1())
    {
        dAi(op1(), n) = dAi(op1(), n) == dAi(op2(), n) ? 1 : 0;
    }
}

void ne(Image<double> &op1, const Image<double> &op2)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op1())
    {
        dAi(op1(), n) = dAi(op1(), n) != dAi(op2(), n) ? 1 : 0;
    }
}

void lt(Image<double> &op1, const Image<double> &op2)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op1())
    {
        dAi(op1(), n) = dAi(op1(), n) < dAi(op2(), n) ? 1 : 0;
    }
}

void le(Image<double> &op1, const Image<double> &op2)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op1())
    {
        dAi(op1(), n) = dAi(op1(), n) <= dAi(op2(), n) ? 1 : 0;
    }
}

void gt(Image<double> &op1, const Image<double> &op2)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op1())
    {
        dAi(op1(), n) = dAi(op1(), n) > dAi(op2(), n) ? 1 : 0;
    }
}

void ge(Image<double> &op1, const Image<double> &op2)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op1())
    {
        dAi(op1(), n) = dAi(op1(), n) >= dAi(op2(), n) ? 1 : 0;
    }
}


double sqrt(Image<double> &op)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op())
    {
        dAi(op(), n) = sqrt(dAi(op(), n));
    }
}

double abs(Image<double> &op)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op())
    {
        dAi(op(), n) = ABS(dAi(op(), n));
    }
}

double log10(Image<double> &op)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op())
    {
        dAi(op(), n) = log10(dAi(op(), n));
    }
}
//4.431-0.4018*LN(ABS(P1+336.6))
double lognikon(Image<double> &op)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op())
    {
        dAi(op(), n) = 4.431-0.4018*log(336.6 + dAi(op(), n));
    }
}

double powerExp = 2;
double power(Image<double> &op)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op())
    {
        dAi(op(), n) = pow(dAi(op(), n), powerExp);
    }
}

int nSlice;
char axis;
double getSlice(Image<double> &op)
{
    MultidimArray<double> imAux;

    op().getSlice(nSlice, imAux, axis);
    op() = imAux;
}


FileName fnOut;
double radialAvg(Image<double> &op)
{
    op().setXmippOrigin();
    Matrix1D<int> center(3);
    center.initZeros();
    MultidimArray<double> radial_mean;
    MultidimArray<int> radial_count;
    radialAverage(op(), center, radial_mean, radial_count);
    radial_mean.write((fnOut.withoutExtension()).addExtension("txt"));

    int my_rad;
    FOR_ALL_ELEMENTS_IN_ARRAY3D(op())
    {
        my_rad = (int)floor(sqrt((double)(i * i + j * j + k * k)));
        op(k, i, j) = radial_mean(my_rad);
    }
}

void ProgOperate::defineParams()
{
    each_image_produces_an_output = true;
    addUsageLine("A simple Xmipp images calculator. Binary and unary operations");
    XmippMetadataProgram::defineParams();
    addParamsLine("== Binary operations: ==");
    addParamsLine("   --plus <file_or_value>    :Sums two images, volumes or adds a numerical value to an image");
    addParamsLine("or --minus <file_or_value>   :Substracts two images, volumes or substracts a numerical value to an image");
    addParamsLine("or --mult <file_or_value>    :Multiplies two images, volumes, or multiplies per a given number");
    addParamsLine("or --divide <file_or_value>  :Divides two images, volumes, or divides per a given number");
    addParamsLine("or --min <file_or_value>     :Minimum of two images, volumes, or number (pixel-wise)");
    addParamsLine("or --max <file_or_value>     :Maximum of two images, volumes, or number (pixel-wise)");
    addParamsLine("or --compare <file_or_value> :Returns -1 if the left value is less, 0 if are equal or 1 if greater.(pixel-wise)");
    addParamsLine("==+ Relational operations: ==");
    addParamsLine("or --eq <file_or_value>      :Returns 1 if the pixels values are equal, 0 otherwise (pixel-wise)");
    addParamsLine("or --le <file_or_value>      :Returns 1 if the pixels values are equal  less, 0 otherwise (pixel-wise)");
    addParamsLine("or --lt <file_or_value>      :Returns 1 if the pixels values are equal  less, 0 otherwise (pixel-wise)");
    addParamsLine("or --ge <file_or_value>      :Returns 1 if the pixels values are equal  less, 0 otherwise (pixel-wise)");
    addParamsLine("or --gt <file_or_value>      :Returns 1 if the pixels values are equal  less, 0 otherwise (pixel-wise)");
    addParamsLine("or --ne <file_or_value>      :Returns 1 if the pixels values are equal  less, 0 otherwise (pixel-wise)");

    addParamsLine("== Unary operations: ==");
    addParamsLine("or --log10                   :Computes the logarithm of an image");
    addParamsLine("or --lognikon                :Computes the logarithm + 1 of an image");
    addParamsLine("or --sqrt                    :Computes the square root of an image");
    addParamsLine("or --abs                     :Computes the absolute value of an image");
    addParamsLine("or --pow <value=2>           :Computes the power of an image");
    addParamsLine("or --slice  <value>          :Extracts a given slice from a volume");
    addParamsLine("or --column <value>          :Extracts a given column from a image or volume");
    addParamsLine("or --row    <value>          :Extracts a given row from a image or volume");
    addParamsLine("or --radial_avg              :Compute the radial average of an image");

    addExampleLine("Sum two volumes and save result", false);
    addExampleLine("xmipp_image_operate -i volume1.vol --plus volume2.vol -o result.vol");
    addExampleLine("Calculate the log10 of an image called example.xmp and store the resulting one in example_log.xmp", false);
    addExampleLine("xmipp_image_operate -i example.xmp --log10 -o example_log.xmp");
    addExampleLine("Calculate the square root of a volume called example.vol and store it in expample_sq.vol", false);
    addExampleLine("xmipp_image_operate -i example.vol --sqrt  -o expample_sq.vol");
    addExampleLine("Extract the slice number 10 of a set of of volumes given in the sel file called volumes.sel. The names of the output images its supposed to be in the selfile images.sel", false);
    addExampleLine("xmipp_image_operate -i volumes.sel --slice 10 -o images.sel");
    addExampleLine("Sum 5 to every image in images.sel and rewrite the input images", false);
    addExampleLine("xmipp_image_operate -i images.sel --plus 5");
    addExampleLine("Substract two volumes:", false);
    addExampleLine("xmipp_image_operate -i volume1.vol --minus volume2.vol -o volume3.vol");
    addExampleLine("Multiply an image by 2 in every pixel:", false);
    addExampleLine("xmipp_image_operate -i image.xmp --mult 2 -o image2.xmp");
    addExampleLine("Divide 2 by the value of every pixel in the image:", false);
    addExampleLine("xmipp_image_operate -i 2 -divide image.xmp -o image2.xmp");
    addExampleLine(" Rotational average", false);
    addExampleLine("xmipp_image_operate -i image.xmp -radial_avg -o image.rad");
    addExampleLine("where image.rad is an ascii file for plotting the radial_averaged profile, image.rad.img a radial_averaged image", false);
}

void ProgOperate::readParams()
{
    XmippMetadataProgram::readParams();
    binaryOperator = NULL;
    unaryOperator = NULL;
    isValue = false;
    // Check operation to do
    //Binary operations
    if (checkParam("--plus"))
    {
        file_or_value = getParam("--plus");
        binaryOperator = plus;
    }
    else if (checkParam("--minus"))
    {
        file_or_value = getParam("--minus");
        binaryOperator = minus;
    }
    else if (checkParam("--mult"))
    {
        file_or_value = getParam("--mult");
        binaryOperator = mult;
    }
    else if (checkParam("--divide"))
    {
        file_or_value = getParam("--divide");
        binaryOperator = divide;
    }
    else if (checkParam("--max"))
    {
        file_or_value = getParam("--max");
        binaryOperator = max;
    }
    else if (checkParam("--min"))
    {
        file_or_value = getParam("--min");
        binaryOperator = min;
    }
    else if (checkParam("--compare"))
    {
        file_or_value = getParam("--compare");
        binaryOperator = compare;
    }
    ///Relational operations
    else if (checkParam("--eq"))
    {
        file_or_value = getParam("--eq");
        binaryOperator = eq;
    }
    else if (checkParam("--ne"))
    {
        file_or_value = getParam("--ne");
        binaryOperator = ne;
    }
    else if (checkParam("--lt"))
    {
        file_or_value = getParam("--lt");
        binaryOperator = lt;
    }
    else if (checkParam("--le"))
    {
        file_or_value = getParam("--le");
        binaryOperator = le;
    }
    else if (checkParam("--gt"))
    {
        file_or_value = getParam("--gt");
        binaryOperator = gt;
    }
    else if (checkParam("--ge"))
    {
        file_or_value = getParam("--ge");
        binaryOperator = ge;
    }
    ///Unary operations
    else if (checkParam("--log10"))
        unaryOperator = log10;
    else if (checkParam("--lognikon"))
        unaryOperator = lognikon;
    else if (checkParam("--sqrt"))
        unaryOperator = sqrt;
    else if (checkParam("--abs"))
        unaryOperator = abs;
    else if (checkParam("--pow"))
    {
        powerExp = getDoubleParam("--pow");
        unaryOperator = power;
    }
    else if (checkParam("--slice"))
    {
        axis = 'Z';
        nSlice = getIntParam("--slice");
        unaryOperator = getSlice;
    }
    else if (checkParam("--column"))
    {
        axis = 'X';
        nSlice = getIntParam("--column");
        unaryOperator = getSlice;
    }
    else if (checkParam("--row"))
    {
        axis = 'Y';
        nSlice = getIntParam("--row");
        unaryOperator = getSlice;
    }
    else if (checkParam("--radial_avg"))
    {
        fnOut = fn_out;
        unaryOperator = radialAvg;
    }
    else if (checkParam("--log10"))
        unaryOperator = log10;
    else if (checkParam("--lognikon"))
        unaryOperator = lognikon;
    else
        REPORT_ERROR(ERR_VALUE_INCORRECT, "No valid operation specified");
    if (binaryOperator != NULL)
    {
        if (!file_or_value.exists())
        {
            isValue = true;
            value = textToFloat(file_or_value);
            img2().resizeNoCopy(zdimOut, ydimOut, xdimOut);
            img2().initConstant(value);
        }
        else
        {
            isValue = false;
            fn2 = file_or_value;
            if (fn2.isMetaData())
            {
                md2.read(fn2);
                if (mdInSize != md2.size())
                    REPORT_ERROR(ERR_MD, "Both metadatas operands should be of same size.");
                md2Iterator = MDIterator(md2);
            }
            else
            {
                isValue = true;
                img2.read(fn2);
            }
        }
    }
}

void ProgOperate::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
    Image<double> img;
    img.readApplyGeo(fnImg, rowIn);

    if (unaryOperator != NULL)
        unaryOperator(img);
    else
    {
        if (!isValue)
        {
            img2.readApplyGeo(md2, md2Iterator.objId);
            md2Iterator.moveNext();
        }
        binaryOperator(img, img2);
    }
    img.write(fnImgOut);
}
