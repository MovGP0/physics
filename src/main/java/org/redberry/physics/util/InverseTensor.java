/*
 *
 * Redberry: symbolic tensor computations library.
 * Copyright (C) 2010-2011  Stanislav Poslavsky <stvlpos@mail.ru>
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
package org.redberry.physics.util;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import redberry.core.tensor.*;
import redberry.core.tensor.generator.GeneratedTensor;
import redberry.core.tensor.generator.TensorGeneratorSP;
import redberry.core.tensor.indexmapping.IndexMappings;
import redberry.core.tensor.test.TTest;
import redberry.core.transformation.CalculateNumbers;
import redberry.core.transformation.Transformation;
import redberry.core.transformation.collect.FullScalarsSplitCriteria;
import redberry.core.transformation.concurrent.ExpandAndCollectTransformation;
import redberry.core.transformation.contractions.IndexesContractionsTransformation;
import redberry.core.utils.Indicator;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class InverseTensor {
    final List<Expression> linearEquations;
    final Tensor[] variables;

    public InverseTensor(Expression toInverse, Expression equation, Tensor[] samples) {
        this(toInverse, equation, samples, new Transformation[0]);
    }

    public InverseTensor(Expression toInverse, Expression equation, Tensor[] samples, Transformation[] transformations) {
        Product leftEq = (Product) equation.left();
        Tensor inverseLhs = null;
        for (Tensor t : leftEq)
            if (!IndexMappings.mappingExists(t, toInverse.left(), true)) {
                inverseLhs = t.clone();
                break;
            }
        GeneratedTensor generatedTensor = TensorGeneratorSP.generateStructure(inverseLhs.getIndexes(), samples);
        variables = generatedTensor.coefficients;
        Expression inverse = new Expression(inverseLhs, generatedTensor.generatedTensor);

        Transformation[] transformations1 = new Transformation[transformations.length + 2];
        System.arraycopy(transformations, 0, transformations1, 0, transformations.length);
        transformations1[transformations.length] =
                IndexesContractionsTransformation.CONTRACTIONS_WITH_METRIC;
        transformations1[transformations.length + 1] = CalculateNumbers.INSTANCE;
        final Transformation expandCollect =
                new ExpandAndCollectTransformation(
                FullScalarsSplitCriteria.INSTANCE,
                Indicator.FALSE_INDICATOR,
                transformations1);

        equation.eval(
                toInverse.asSubstitution(),
                inverse.asSubstitution(),
                expandCollect);


        Split rightSplit = split(equation.right());
        linearEquations = new ArrayList<>();
        for (Tensor summand : equation.left()) {
            Split current = split(summand);
            if (TTest.testParity(current.nonScalar, rightSplit.nonScalar))
                linearEquations.add(new Expression(current.scalar, rightSplit.scalar));
            else
                linearEquations.add(new Expression(current.scalar, TensorNumber.createZERO()));
        }
        generateMapleFile();

    }

    private static class Split {
        final Tensor scalar;
        final Tensor nonScalar;

        public Split(Tensor scalar, Tensor nonScalar) {
            this.scalar = scalar;
            this.nonScalar = nonScalar;
        }
    }

    private static Split split(Tensor tensor) {
        Tensor scalar = null, nonScalar = null;
        if (tensor instanceof Product) {
            ProductContent pc = (ProductContent) tensor.getContent();
            if (pc.getScalarContents().length == 0) {
                nonScalar = tensor;
                scalar = TensorNumber.createONE();
            } else {
                scalar = new Product();
                if (!pc.getFactor().isOne())
                    ((Product) scalar).add(new TensorNumber(pc.getFactor()));
                for (ProductContent content : pc.getScalarContents())
                    ((Product) scalar).add(content.getRange(0, content.size()));
                nonScalar = new Product(pc.getNonScalarContent().
                        getRange(0, pc.getNonScalarContent().size()));
            }
        } else {
            nonScalar = tensor;
            scalar = TensorNumber.createONE();
        }
        return new Split(scalar, nonScalar);
    }

    public void generateMapleFile() {
        try {
            FileOutputStream output = new FileOutputStream("equations.maple");
            PrintStream file = new PrintStream(output);

            file.append("ans:=array([");
            for (int i = 0; i < variables.length; ++i)
                if (i == variables.length - 1)
                    file.append(variables[i].toString());
                else
                    file.append(variables[i] + ",");
            file.append("]);\n");

            file.println("eq:=array(1.." + linearEquations.size() + ");");
            for (int i = 0; i < linearEquations.size(); i++)
                file.println("eq[" + (i + 1) + "]:=" + linearEquations.get(i) + ";");

            //for (int i = 0; i < systemEquations.size() - 1; i++)
            //   file.append("eq" + (i + 1) + ",");
            file.print("Result := solve({seq(eq[i],i=1.." + linearEquations.size() + ")},[");
            //file.append("eq" + systemEquations.size() + "},[");
            //file.append(",[");
            for (int i = 0; i < variables.length; ++i)
                if (i == variables.length - 1)
                    file.append(variables[i].toString());
                else
                    file.append(variables[i] + ",");
            file.append("]);\n");
            
            file.println("file:=fopen(\"/home/stas/NetBeansProjects/redberry-physics/" + "equations.mapleOut\",WRITE);");
            //    file.append("try\n");
            file.append("fprintf(file,\"ans" + variables.length + ":=array(1.." + variables.length + ");\\n\");\n");
            file.append("for k from 1 to " + variables.length + " do\n");
            file.append("fprintf(file,\"ans" + variables.length + "[%a]:=%a=%a;\\n\",k,ans" + variables.length + "[k],rhs(Result[1][k]))\n");
            file.append("od;\n");
            // file.append("finally fclose(file) end try;");
            file.append("fclose(file);");
        } catch (Exception e) {
        }
    }
}
