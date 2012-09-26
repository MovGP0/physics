/*
 * Redberry: symbolic tensor computations.
 *
 * Copyright (c) 2010-2012:
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

import cc.redberry.core.indexmapping.IndexMappings;
import cc.redberry.core.number.*;
import cc.redberry.core.tensor.*;
import cc.redberry.core.tensor.iterator.*;
import cc.redberry.core.tensorgenerator.*;
import cc.redberry.core.transformations.*;
import cc.redberry.core.utils.*;
import java.io.*;
import java.util.*;
import org.apache.commons.math3.random.*;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public final class InverseTensor {

    private final Expression[] equations;
    private final SimpleTensor[] uncknownCoefficients;
    private final Expression generalInverse;

    public InverseTensor(Expression toInverse, Expression equation, Tensor[] samples) {
        this(toInverse, equation, null, false, new Transformation[0]);
    }

    public InverseTensor(Expression toInverse, Expression equation, Tensor[] samples, boolean symmetricForm, Transformation[] transformations) {
        Product leftEq = (Product) equation.get(0);
        Tensor inverseLhs = null;
        for (Tensor t : leftEq)
            if (!IndexMappings.mappingExists(t, toInverse.get(0))) {
                inverseLhs = t;
                break;
            }
        GeneratedTensor generatedTensor = TensorGenerator.generateStructure("c", inverseLhs.getIndices(), symmetricForm, samples);
        uncknownCoefficients = generatedTensor.coefficients;
        generalInverse = ExpressionFactory.FACTORY.create(inverseLhs, generatedTensor.generatedTensor);

        Tensor temp = equation;
        temp = toInverse.transform(temp);
        temp = generalInverse.transform(temp);
        transformations = ArraysUtils.addAll(new Transformation[]{ContractIndices.ContractIndices}, transformations);
        temp = Expand.expand(temp, transformations);

        for (Transformation transformation : transformations)
            temp = transformation.transform(temp);
        temp = CollectNonScalars.collectNonScalars(temp);
        equation = (Expression) temp;

        List<Split> rightSplit = new ArrayList<>();

        if (equation.get(1) instanceof Sum)
            for (Tensor summand : equation.get(1))
                rightSplit.add(Split.splitScalars(summand));
        else
            rightSplit.add(Split.splitScalars(equation.get(1)));
        List<Expression> equationsList = new ArrayList<>();
        for (Tensor summand : equation.get(0)) {
            Split current = Split.splitScalars(summand);
            boolean one = false;
            for (Split split : rightSplit)
                if (TensorUtils.equals(current.factor, split.factor)) {
                    equationsList.add(ExpressionFactory.FACTORY.create(current.summand, split.summand));
                    one = true;
                    break;
                }
            if (!one)
                equationsList.add(ExpressionFactory.FACTORY.create(current.summand, Complex.ZERO));
        }
        this.equations = equationsList.toArray(new Expression[equationsList.size()]);
    }

    public Expression[] getEquations() {
        return equations.clone();
    }

    public Expression getGeneralInverseForm() {
        return generalInverse;
    }

    public SimpleTensor[] getUncknownCoefficients() {
        return uncknownCoefficients.clone();
    }

    public static Tensor findInverseWithMaple(Expression toInverse,
                                              Expression equation,
                                              Tensor[] samples,
                                              boolean symmetricForm,
                                              Transformation[] transformations,
                                              String mapleBinDir,
                                              String path) {
        return findInverseWithMaple(toInverse, equation, samples, symmetricForm, false, transformations, mapleBinDir, path);
    }

    public static Tensor findInverseWithMaple(Expression toInverse,
                                              Expression equation,
                                              Tensor[] samples,
                                              boolean symmetricForm,
                                              boolean keepFreeParametres,
                                              Transformation[] transformations,
                                              String mapleBinDir,
                                              String path) {
        InverseTensor inverseTensor = new InverseTensor(toInverse, equation, samples, symmetricForm, transformations);
        final Expression[] equations = inverseTensor.equations.clone();

        HashMap<TensorWrapperWithEquals, Tensor> tensorSubstitutions = new HashMap<>();
        SymbolsGenerator generator = new SymbolsGenerator("scalar", ArraysUtils.addAll(samples, toInverse, equation));
        int i;
        for (i = 0; i < equations.length; ++i) {
            Expression eq = equations[i];
            TensorLastIterator iterator = new TensorLastIterator(eq);
            Tensor t;
            while ((t = iterator.next()) != null) {
                if (!(t instanceof Product) || t.getIndices().size() == 0)
                    continue;
                Split split = Split.splitIndexless(t);
                if (split.factor.getIndices().size() == 0)
                    continue;
                TensorWrapperWithEquals twwe = new TensorWrapperWithEquals(split.factor);
                if (!tensorSubstitutions.containsKey(twwe)) {
                    Tensor s;
                    tensorSubstitutions.put(twwe, s = generator.take());
                    iterator.set(Tensors.multiply(s, split.summand));
                } else
                    iterator.set(Tensors.multiply(tensorSubstitutions.get(twwe), split.summand));
            }
            equations[i] = (Expression) iterator.result();
        }

        System.out.println("Inverse tensor: " + inverseTensor.generalInverse);
        System.out.println();
        try {
            FileOutputStream output = new FileOutputStream(path + "/equations.maple");
            PrintStream file = new PrintStream(output);
            file.append("with(StringTools):\n");
            file.append("ans:=array([");
            for (i = 0; i < inverseTensor.uncknownCoefficients.length; ++i)
                if (i == inverseTensor.uncknownCoefficients.length - 1)
                    file.append(inverseTensor.uncknownCoefficients[i].toString());
                else
                    file.append(inverseTensor.uncknownCoefficients[i] + ",");
            file.append("]):\n");

            file.println("eq:=array(1.." + equations.length + "):");
            for (i = 0; i < equations.length; i++)
                file.println("eq[" + (i + 1) + "]:=" + equations[i] + ":");

            file.print("Result := solve({seq(eq[i],i=1.." + equations.length + ")},[");
            for (i = 0; i < inverseTensor.uncknownCoefficients.length; ++i)
                if (i == inverseTensor.uncknownCoefficients.length - 1)
                    file.append(inverseTensor.uncknownCoefficients[i].toString());
                else
                    file.append(inverseTensor.uncknownCoefficients[i] + ",");
            file.append("]);\n");
            file.println("file:=fopen(\"" + path + "/equations.mapleOut\",WRITE):");
            file.append("if nops(Result) <> 0 then\n");
            file.append("for k from 1 to " + inverseTensor.uncknownCoefficients.length + " do\n");
            file.append("temp1 := SubstituteAll(convert(lhs(Result[1][k]), string), \"^\", \"**\");\n");
            file.append("temp2 := SubstituteAll(convert(rhs(Result[1][k]), string), \"^\", \"**\");\n");
            file.append("fprintf(file,\"%s=%s\\n\",temp1,temp2);\n");
            file.append("od:\n");
            file.append("end if;\n");
            file.append("fclose(file):");
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
        
        //TODO use temp file
        //File.createTempFile(path, path)
        
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

        try {
            Expression[] coefficientsResults = new Expression[inverseTensor.uncknownCoefficients.length];
            FileInputStream fstream = new FileInputStream(path + "/equations.mapleOut");
            if (fstream.available() == 0)
                return null;
            DataInputStream in = new DataInputStream(fstream);
            BufferedReader br = new BufferedReader(new InputStreamReader(in));
            String strLine;
            i = -1;
            while ((strLine = br.readLine()) != null)
                coefficientsResults[++i] = Tensors.parseExpression(strLine);
            Tensor inverse = inverseTensor.generalInverse;
            for (Expression coef : coefficientsResults)
                if (coef.isIdentity()) {
                    if (!keepFreeParametres)
                        inverse = Tensors.expression(coef.get(0), Complex.ZERO).transform(inverse);
                } else
                    inverse = (Expression) coef.transform(inverse);
            for (Map.Entry<TensorWrapperWithEquals, Tensor> entry : tensorSubstitutions.entrySet())
                inverse = Tensors.expression(entry.getValue(), entry.getKey().tensor).transform(inverse);
            in.close();
            return inverse;
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    private static final class TensorWrapperWithEquals {

        private final Tensor tensor;

        public TensorWrapperWithEquals(Tensor tensor) {
            this.tensor = tensor;
        }

        @Override
        public boolean equals(Object obj) {
            if (obj == null)
                return false;
            if (getClass() != obj.getClass())
                return false;
            final TensorWrapperWithEquals other = (TensorWrapperWithEquals) obj;
            return TensorUtils.equals(tensor, other.tensor);
        }

        @Override
        public int hashCode() {
            return tensor.hashCode();
        }
    }
}