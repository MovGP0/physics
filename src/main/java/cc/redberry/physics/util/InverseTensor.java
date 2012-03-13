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
package cc.redberry.physics.util;

import cc.redberry.core.combinatorics.Symmetries;
import cc.redberry.core.context.CC;
import cc.redberry.core.indexmapping.IndexMappings;
import cc.redberry.core.tensor.*;
import cc.redberry.core.tensor.testing.TTest;
import cc.redberry.core.utils.Indicator;
import cc.redberry.tensorgenerator.GeneratedTensor;
import cc.redberry.tensorgenerator.TensorGenerator;
import cc.redberry.transformation.CalculateNumbers;
import cc.redberry.transformation.Transformation;
import cc.redberry.transformation.collect.FullScalarsSplitCriteria;
import cc.redberry.transformation.concurrent.ExpandAndCollectTransformation;
import cc.redberry.transformation.contractions.IndicesContractionsTransformation;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;


/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class InverseTensor {
    final List<Expression> linearEquations;
    final Tensor[] variables;
    final Expression inverse;

    public InverseTensor(Expression toInverse, Expression equation, Tensor[] samples) {
        this(toInverse, equation, null, samples, new Transformation[0]);
    }

    public InverseTensor(Expression toInverse, Expression equation, Tensor[] samples, Transformation[] transformations) {
        this(toInverse, equation, null, samples, transformations);
    }

    public InverseTensor(Expression toInverse, Expression equation, Symmetries symmeties, Tensor[] samples, Transformation[] transformations) {
        Product leftEq = (Product) equation.left();
        Tensor inverseLhs = null;
        for (Tensor t : leftEq)
            if (!IndexMappings.mappingExists(t, toInverse.left(), true)) {
                inverseLhs = t.clone();
                break;
            }
        GeneratedTensor generatedTensor = TensorGenerator.generateStructure(inverseLhs.getIndices(), symmeties, samples);
        variables = generatedTensor.coefficients;
        inverse = new Expression(inverseLhs, generatedTensor.generatedTensor);

        Transformation[] transformations1 = new Transformation[transformations.length + 2];
        transformations1[0] =
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC;
        System.arraycopy(transformations, 0, transformations1, 1, transformations.length);
        transformations1[transformations.length + 1] = CalculateNumbers.INSTANCE;
        final Transformation expandCollect =
                new ExpandAndCollectTransformation(
                FullScalarsSplitCriteria.INSTANCE,
                Indicator.SYMBOL_INDICATOR,
                transformations1);

        equation.eval(
                toInverse.asSubstitution(),
                inverse.asSubstitution(),
                expandCollect);

        List<Split> rightSplit = new ArrayList<>();

        if (equation.right() instanceof Sum)
            for (Tensor summand : equation.right())
                rightSplit.add(split(summand));
        else
            rightSplit.add(split(equation.right()));
        linearEquations = new ArrayList<>();
        for (Tensor summand : equation.left()) {
            Split current = split(summand);
            boolean one = false;
            for (Split split : rightSplit)
                if (TTest.testParity(current.nonScalar, split.nonScalar)) {
                    linearEquations.add(new Expression(current.scalar, split.scalar));
                    one = true;
                    break;
                }
            if (!one)
                linearEquations.add(new Expression(current.scalar, TensorNumber.createZERO()));
        }
        if (CC.getMapleDirectory() != null)
            generateMapleFile();
    }

    public Expression inverse() {
        return inverse;
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
                nonScalar = new Product(pc.getRange(0, pc.size()));
                scalar = new TensorNumber(pc.getFactor());
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
        return new Split(scalar.equivalent(), nonScalar.equivalent());
    }

    private void generateMapleFile() {
        String workingFolder = CC.getWorkingFolder() + "/" + InverseTensor.class.getSimpleName();
        try {
            (new File(workingFolder)).mkdir();
        } catch (Exception e) {
            System.err.println("Error: cannot create working folder" + e.getMessage());
        }

        String mapleDirectory = CC.getMapleDirectory();
        try {
            FileOutputStream output = new FileOutputStream(workingFolder + "/equations.maple");
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

            file.print("Result := solve({seq(eq[i],i=1.." + linearEquations.size() + ")},[");
            for (int i = 0; i < variables.length; ++i)
                if (i == variables.length - 1)
                    file.append(variables[i].toString());
                else
                    file.append(variables[i] + ",");
            file.append("]);\n");

            file.println("file:=fopen(\"" + workingFolder + "/equations.mapleOut\",WRITE);");
            file.append("for k from 1 to " + variables.length + " do\n");
            file.append("fprintf(file,\"%a=%a;\\n\",lhs(Result[1][k]),rhs(Result[1][k]))\n");
            file.append("od;\n");
            file.append("fclose(file);");
        } catch (Exception e) {
        }
    }
}
