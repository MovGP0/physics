/*
 * Redberry: symbolic tensor computations.
 *
 * Copyright (c) 2010-2013:
 *   Stanislav Poslavsky   <stvlpos@mail.ru>
 *   Bolotin Dmitriy       <bolotin.dmitriy@gmail.com>
 *
 * This file is part of Redberry.
 *
 * Redberry is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Redberry is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Redberry. If not, see <http://www.gnu.org/licenses/>.
 */

package cc.redberry.physics.utils;

import cc.redberry.core.context.CC;
import cc.redberry.core.indexmapping.IndexMappings;
import cc.redberry.core.number.Complex;
import cc.redberry.core.tensor.*;
import cc.redberry.core.tensor.iterator.FromChildToParentIterator;
import cc.redberry.core.tensorgenerator.GeneratedTensor;
import cc.redberry.core.tensorgenerator.SymbolsGenerator;
import cc.redberry.core.tensorgenerator.TensorGenerator;
import cc.redberry.core.transformations.CollectNonScalarsTransformation;
import cc.redberry.core.transformations.EliminateMetricsTransformation;
import cc.redberry.core.transformations.expand.ExpandTransformation;
import cc.redberry.core.transformations.Transformation;
import cc.redberry.core.utils.ArraysUtils;
import cc.redberry.core.utils.THashMap;
import cc.redberry.core.utils.TensorUtils;

import java.io.*;
import java.util.*;

/**
 * This class provides opportunities to find inverse of tensor. In
 * other words it can be used to solve the the equation of the form
 * <pre>
 *     <i>T^{ij..}_{kp..}*Tinv^{kp..}_{mn..} = d^{i}_{m}*d^{j}_{m}*.. + ... (combinations of kroneckers),</i>
 * </pre>
 * where <i>T</i> is specified tensor, <i>Tinv</i> is unknown tensor
 * and <i>d</i> - kronecker delta.
 * <p/>
 * The main goal of this class is to create the tensor <i>Tinv</i> of the
 * most general form with unknown coefficients, and produce a system of
 * linear equations on these coefficients. The resulting equations can
 * then be solved and coefficients values substituted in generated
 * <i>Tinv</i> tensor.
 * </p>
 * <p/>
 * <br>The following example demonstrates the usage of the {@code InverseTensor} to
 * find out the photon propagator in Lorentz gauge:</br>
 * <pre>
 *      ...
 *      //expression specifies tensor, which need to inverse
 *      Expression toInverse = Tensors.parseExpression("D_mn = k_m*k_n-(1/a)*k_i*k^i*g_mn");
 *      //linear equation on the unknown tensor K
 *      Expression equation = Tensors.parseExpression("D_ab*K^ac=d_b^c");
 *      //samples from which inverse should be formed
 *      Tensor[] samples = {Tensors.parse("g_mn"), Tensors.parse("g^mn"), Tensors.parse("d_m^n"), Tensors.parse("k_m"), Tensors.parse("k^b")};
 *
 *      InverseTensor inverseTensor = new InverseTensor(toInverse,equation,samples);
 *      System.out.println(inverseTensor.getGeneralInverseForm());
 *      System.out.println(Arrays.toString(inverseTensor.getEquations()));
 * </pre>
 * <br>The above code displays the inverse of specified tensor</br>
 * <pre>
 *     K^{ac} = a1*g^{ac}+a0*k^{a}*k^{c}
 * </pre>
 * and a system of equations on its coefficients
 * <pre>
 *     [(-a**(-1)*a0+a0)*k^{i}*k_{i}+a1 = 0, -a**(-1)*a1*k_{i}*k^{i} = 1]
 * </pre>
 * </p>
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 * @see #findInverseWithMaple(cc.redberry.core.tensor.Expression, cc.redberry.core.tensor.Expression, cc.redberry.core.tensor.Tensor[], boolean, boolean, cc.redberry.core.transformations.Transformation[], String, String)
 * @since 1.0
 */
public final class InverseTensor {

    private final Expression[] equations;
    private final SimpleTensor[] unknownCoefficients;
    private final Expression generalInverse;

    /**
     * Creates the {@code InverseTensor} instance from the equation.
     *
     * @param toInverse expression specifies tensor, which need to inverse
     * @param equation  linear equation on the unknown tensor in the form
     *                  T^{..}_{...}*Tinv^{...}_{...} = ...
     * @param samples   samples from which inverse should be formed
     *                  into account when forming a system of linear equations
     */
    public InverseTensor(Expression toInverse, Expression equation, Tensor[] samples) {
        this(toInverse, equation, samples, false, new Transformation[0]);
    }

    /**
     * Creates the {@code InverseTensor} instance from the equation.
     *
     * @param toInverse       expression specifies tensor, which need to inverse
     * @param equation        linear equation on the unknown tensor in the form
     *                        T^{..}_{...}*Tinv^{...}_{...} = ...
     * @param samples         samples from which inverse should be formed
     * @param symmetricForm   specifies whether inverse tensor should be symmetric
     * @param transformations additional simplification rules, which can be taken
     *                        into account when forming a system of linear equations
     */
    public InverseTensor(Expression toInverse,
                         Expression equation,
                         Tensor[] samples,
                         boolean symmetricForm,
                         Transformation[] transformations) {
        if (!(equation.get(0) instanceof Product))
            throw new IllegalArgumentException("Equation l.h.s. is not a product of tensors.");

        Product leftEq = (Product) equation.get(0);

        //matching toInverse l.h.s in equation
        Tensor inverseLhs = null;
        for (Tensor t : leftEq)
            if (!IndexMappings.mappingExists(t, toInverse.get(0))) {
                inverseLhs = t;
                break;
            }

        //creating tensor of the most general form from the specified samples
        GeneratedTensor generatedTensor = TensorGenerator.generateStructure(newCoefficientName(toInverse, equation), inverseLhs.getIndices(), symmetricForm, samples);
        unknownCoefficients = generatedTensor.coefficients;
        //creating inverse tensor expression
        generalInverse = Tensors.expression(inverseLhs, generatedTensor.generatedTensor);

        //substituting toInverse and generalInverse into equation
        Tensor temp = equation;
        temp = toInverse.transform(temp);
        temp = generalInverse.transform(temp);

        //collecting all transformations in single array
        transformations = ArraysUtils.addAll(new Transformation[]{EliminateMetricsTransformation.ELIMINATE_METRICS}, transformations);

        //preparing equation
        temp = ExpandTransformation.expand(temp, transformations);
        for (Transformation transformation : transformations)
            temp = transformation.transform(temp);
        temp = CollectNonScalarsTransformation.collectNonScalars(temp);
        equation = (Expression) temp;


        //processing r.h.s. of the equation
        List<Split> rightSplit = new ArrayList<>();
        if (equation.get(1) instanceof Sum)
            for (Tensor summand : equation.get(1))
                rightSplit.add(Split.splitScalars(summand));
        else
            rightSplit.add(Split.splitScalars(equation.get(1)));

        //forming system of linear equations
        List<Expression> equationsList = new ArrayList<>();
        for (Tensor summand : equation.get(0)) {
            Split current = Split.splitScalars(summand);
            boolean one = false;
            for (Split split : rightSplit)
                if (TensorUtils.equals(current.factor, split.factor)) {
                    equationsList.add(Tensors.expression(current.summand, split.summand));
                    one = true;
                    break;
                }
            if (!one)
                equationsList.add(Tensors.expression(current.summand, Complex.ZERO));
        }
        this.equations = equationsList.toArray(new Expression[equationsList.size()]);
    }

    private static String newCoefficientName(Tensor... tensors) {
        Set<SimpleTensor> simpleTensors = TensorUtils.getAllSymbols(tensors);
        List<Character> forbidden = new ArrayList<>();
        for (SimpleTensor tensor : simpleTensors) {
            String name = CC.getNameDescriptor(tensor.getName()).getName(tensor.getIndices());
            try {
                Integer.parseInt(name.substring(1));
                forbidden.add(name.charAt(0));
            } catch (NumberFormatException e) {
            }
        }
        Collections.sort(forbidden);
        char c = 'a';
        for (int i = 0; i < forbidden.size(); ++i) {
            if (c != forbidden.get(i).charValue())
                break;
            else {
                ++c;
            }
        }
        return String.valueOf(c);
    }

    /**
     * Return the resulting equations on the unknown coefficients.
     *
     * @return the resulting equations on the unknown coefficients
     */
    public Expression[] getEquations() {
        return equations.clone();
    }

    /**
     * Returns the inverse of the tensor with unknown coefficients.
     *
     * @return the inverse of the tensor with unknown coefficients
     */
    public Expression getGeneralInverseForm() {
        return generalInverse;
    }

    /**
     * Returns the array of the unknown coefficients.
     *
     * @return the array of the unknown coefficients
     */
    public SimpleTensor[] getUnknownCoefficients() {
        return unknownCoefficients.clone();
    }

    /**
     * This method calculates the tensor inverse to the specified tensor according
     * to the specified equation using the Maple facilities to solve the system of
     * linear equations.  The Maple code will be placed in the specified temporary
     * directory in {@code equations.maple} file. The solution of the linear system,
     * produced by Maple, will be placed in the specified temporary directory in
     * {@code equations.mapleOut} file.
     * <p/>
     * <br>The following example demonstrates the usage of this method to
     * find out the photon propagator in Lorentz gauge:</br>
     * <pre>
     *      ...
     *      //expression specifies tensor, which need to inverse
     *      Expression toInverse = Tensors.parseExpression("D_mn = k_m*k_n-(1/a)*k_i*k^i*g_mn");
     *      //linear equation on the unknown tensor K
     *      Expression equation = Tensors.parseExpression("D_ab*K^ac=d_b^c");
     *      //samples from which inverse should be formed
     *      Tensor[] samples = {Tensors.parse("g_mn"), Tensors.parse("g^mn"), Tensors.parse("d_m^n"), Tensors.parse("k_m"), Tensors.parse("k^b")};
     *
     *      Tensor inverse = InverseTensor.findInverseWithMaple(toInverse, equation, samples, false, new Transformation[0], mapleBinDir, temporaryDir);
     *      System.out.println(inverse);
     * </pre>
     * <br>The above code displays the inverse of specified tensor</br>
     * <pre>
     *     K^ac=-a*g^ac*(k_i*k^i)**(-1)+a**2/(a-1)*k^a*k^c*(k_i*k^i)**(-2)
     * </pre>
     * <p/>
     *
     * @param toInverse       expression specifies tensor, which need to inverse
     * @param equation        linear equation on the unknown tensor in the form
     *                        T^{..}_{...}*Tinv^{...}_{...} = ...
     * @param samples         samples from which inverse should be formed
     * @param symmetricForm   specifies whether inverse tensor should be symmetric
     * @param transformations additional simplification rules, which can be taken
     *                        into account when forming a system of linear equations
     * @param mapleBinDir     path to Maple bin directory (e.g. "/home/user/maple14/bin")
     * @param path            path to your temporary folder
     * @return tensor inverse to the specified tensor according to the specified
     *         equation and null if inverse does not exist
     */
    public static Tensor findInverseWithMaple(Expression toInverse,
                                              Expression equation,
                                              Tensor[] samples,
                                              boolean symmetricForm,
                                              Transformation[] transformations,
                                              String mapleBinDir,
                                              String path) {
        return findInverseWithMaple(toInverse, equation, samples, symmetricForm, false, transformations, mapleBinDir, path);
    }

    /**
     * This method calculates the tensor inverse to the specified tensor according
     * to the specified equation using the Maple facilities to solve the system of
     * linear equations. The Maple code will be placed in the specified temporary
     * directory in {@code equations.maple} file. The solution of the linear system,
     * produced by Maple, will be placed in the specified temporary directory in
     * {@code equations.mapleOut} file.
     * <p/>
     * <p/>
     * <br>The following example demonstrates the usage of this method to
     * find out the photon propagator in Lorentz gauge:</br>
     * <pre>
     *      ...
     *      //expression specifies tensor, which need to inverse
     *      Expression toInverse = Tensors.parseExpression("D_mn = k_m*k_n-(1/a)*k_i*k^i*g_mn");
     *      //linear equation on the unknown tensor K
     *      Expression equation = Tensors.parseExpression("D_ab*K^ac=d_b^c");
     *      //samples from which inverse should be formed
     *      Tensor[] samples = {Tensors.parse("g_mn"), Tensors.parse("g^mn"), Tensors.parse("d_m^n"), Tensors.parse("k_m"), Tensors.parse("k^b")};
     *
     *      Tensor inverse = InverseTensor.findInverseWithMaple(toInverse, equation, samples, false, new Transformation[0], false, mapleBinDir, temporaryDir);
     *      System.out.println(inverse);
     * </pre>
     * <br>The above code displays the inverse of specified tensor</br>
     * <pre>
     *     K^ac=-a*g^ac*(k_i*k^i)**(-1)+a**2/(a-1)*k^a*k^c*(k_i*k^i)**(-2)
     * </pre>
     * <p/>
     * If the produced  system of linear equations have infinitely many solutions, some
     * of the coefficient cannot be determined exactly and remain as free parameters. The
     * flag {@code keepFreeParameters} specifies whether this free parameters should be
     * zeroed.
     *
     * @param toInverse          expression specifies tensor, which need to inverse
     * @param equation           linear equation on the unknown tensor in the form
     *                           T^{..}_{...}*Tinv^{...}_{...} = ...
     * @param samples            samples from which inverse should be formed
     * @param symmetricForm      specifies whether inverse tensor should be symmetric
     * @param transformations    additional simplification rules, which can be taken
     *                           into account when forming a system of linear equations
     * @param keepFreeParameters specifies whether the free parameters remaining from solution
     *                           of linear system should be zeroed
     * @param mapleBinDir        path to Maple bin directory (e.g. "/home/user/maple14/bin")
     * @param path               path to your temporary folder
     * @return tensor inverse to the specified tensor according to the specified
     *         equation and null if inverse does not exist
     */
    public static Tensor findInverseWithMaple(Expression toInverse,
                                              Expression equation,
                                              Tensor[] samples,
                                              boolean symmetricForm,
                                              boolean keepFreeParameters,
                                              Transformation[] transformations,
                                              String mapleBinDir,
                                              String path) {
        //create the general form of the inverse and system of linear equations
        InverseTensor inverseTensor = new InverseTensor(toInverse, equation, samples, symmetricForm, transformations);
        final Expression[] equations = inverseTensor.equations.clone();

        /*in order to process equations with Maple we must to replace all tensors
        with indices (they are found only in scalar combinations) with some symbols*/

        //scalar tensor <-> symbol
        THashMap<Tensor, Tensor> tensorSubstitutions = new THashMap<>();
        //all symbols will have names scalar1,scalar2, etc.
        SymbolsGenerator generator = new SymbolsGenerator("scalar", ArraysUtils.addAll(samples, toInverse, equation));

        //processing equations
        int i;
        for (i = 0; i < equations.length; ++i) {
            Expression eq = equations[i];
            //iterating over the whole equation
            FromChildToParentIterator iterator = new FromChildToParentIterator(eq);
            Tensor t;
            while ((t = iterator.next()) != null) {
                if (!(t instanceof Product) || t.getIndices().size() == 0)
                    continue;
                //split tensor into symbolic and tensor (with nonzero length of indices) parts
                Split split = Split.splitIndexless(t);
                if (split.factor.getIndices().size() == 0)//there is no nonsymbolic part in tensor, e.g. t = a*b
                    continue;

                //there is non symbolic part in current tensor, e.g. t = a*k_{i}*k^{i}, so
                //split.summand = a, split.factor = k_{i}*k^{i}

                if (!tensorSubstitutions.containsKey(split.factor)) {
                    //map does not contains rule for current scalar (e.g. k_{i}*k^{i})
                    Tensor s;
                    //adding new rule for the scalar, e.g. k_{i}*k^{i} = scalar2
                    tensorSubstitutions.put(split.factor, s = generator.take());
                    //replacing this scalar with symbol
                    iterator.set(Tensors.multiply(s, split.summand));
                } else
                    //map is already contains rule for current scalar
                    //replacing this scalar with symbol from map
                    iterator.set(Tensors.multiply(tensorSubstitutions.get(split.factor), split.summand));
            }
            equations[i] = (Expression) iterator.result();
        }

        System.out.println("Inverse tensor: " + inverseTensor.generalInverse);
        System.out.println();

        //creating file with Maple code to solve the system of equations
        try {
            FileOutputStream output = new FileOutputStream(path + "/equations.maple");
            PrintStream file = new PrintStream(output);
            file.append("with(StringTools):\n");
            file.append("ans:=array([");
            for (i = 0; i < inverseTensor.unknownCoefficients.length; ++i)
                if (i == inverseTensor.unknownCoefficients.length - 1)
                    file.append(inverseTensor.unknownCoefficients[i].toString());
                else
                    file.append(inverseTensor.unknownCoefficients[i] + ",");
            file.append("]):\n");

            file.println("eq:=array(1.." + equations.length + "):");
            for (i = 0; i < equations.length; i++)
                file.println("eq[" + (i + 1) + "]:=" + equations[i] + ":");

            file.print("Result := solve({seq(eq[i],i=1.." + equations.length + ")},[");
            for (i = 0; i < inverseTensor.unknownCoefficients.length; ++i)
                if (i == inverseTensor.unknownCoefficients.length - 1)
                    file.append(inverseTensor.unknownCoefficients[i].toString());
                else
                    file.append(inverseTensor.unknownCoefficients[i] + ",");
            file.append("]);\n");
            file.println("file:=fopen(\"" + path + "/equations.mapleOut\",WRITE):");
            file.append("if nops(Result) <> 0 then\n");
            file.append("for k from 1 to " + inverseTensor.unknownCoefficients.length + " do\n");
            file.append("temp1 := SubstituteAll(convert(lhs(Result[1][k]), string), \"^\", \"**\");\n");
            file.append("temp2 := SubstituteAll(convert(rhs(Result[1][k]), string), \"^\", \"**\");\n");
            file.append("fprintf(file,\"%s=%s\\n\",temp1,temp2);\n");
            file.append("od:\n");
            file.append("end if;\n");
            file.append("fclose(file):");
        } catch (Exception e) {
            throw new RuntimeException(e);
        }

        //TODO maybe use temp file

        //running Maple
        try {
            Process p = Runtime.getRuntime().exec(mapleBinDir + "/maple " + path + "/equations.maple");
            BufferedReader bri = new BufferedReader(new InputStreamReader(p.getInputStream()));
            BufferedReader bre = new BufferedReader(new InputStreamReader(p.getErrorStream()));
            String line;
            while ((line = bri.readLine()) != null)
                System.out.println(line);
            bri.close();
            while ((line = bre.readLine()) != null)
                System.out.println(line);
            bre.close();
            p.waitFor();
        } catch (IOException | InterruptedException ex) {
            throw new RuntimeException(ex);
        }

        //reading the Maple output with the solution
        try {
            //allocating resulting coefficients array
            Expression[] coefficientsResults = new Expression[inverseTensor.unknownCoefficients.length];
            FileInputStream fstream = new FileInputStream(path + "/equations.mapleOut");
            if (fstream.available() == 0)
                return null;
            DataInputStream in = new DataInputStream(fstream);
            BufferedReader br = new BufferedReader(new InputStreamReader(in));
            String strLine;
            i = -1;
            //reading resulting solutions from file
            while ((strLine = br.readLine()) != null)
                coefficientsResults[++i] = Tensors.parseExpression(strLine);
            Tensor inverse = inverseTensor.generalInverse;

            //substituting coefficients into general inverse form
            for (Expression coef : coefficientsResults)
                if (coef.isIdentity())//if current coefficient is free parametr
                {
                    if (!keepFreeParameters)
                        inverse = Tensors.expression(coef.get(0), Complex.ZERO).transform(inverse);
                } else
                    inverse = (Expression) coef.transform(inverse);

            //substituting the renamed tensors combinations
            for (Map.Entry<Tensor, Tensor> entry : tensorSubstitutions.entrySet())
                inverse = Tensors.expression(entry.getValue(), entry.getKey()).transform(inverse);
            in.close();
            return inverse;
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }
}
