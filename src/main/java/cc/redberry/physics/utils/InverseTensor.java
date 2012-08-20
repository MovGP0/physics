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

import cc.redberry.core.combinatorics.symmetries.Symmetries;
import cc.redberry.core.context.CC;
import cc.redberry.core.indexmapping.IndexMappings;
import cc.redberry.core.number.*;
import cc.redberry.core.tensor.*;
import cc.redberry.core.tensorgenerator.GeneratedTensor;
import cc.redberry.core.tensorgenerator.TensorGenerator;
import cc.redberry.core.transformations.*;
import cc.redberry.core.utils.*;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public final class InverseTensor {

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
        Product leftEq = (Product) equation.get(0);
        Tensor inverseLhs = null;
        for (Tensor t : leftEq)
            if (!IndexMappings.mappingExists(t, toInverse.get(0))) {
                inverseLhs = t;
                break;
            }
        GeneratedTensor generatedTensor = TensorGenerator.generateStructure("c", inverseLhs.getIndices(), false, samples);
        variables = generatedTensor.coefficients;
        inverse = ExpressionFactory.FACTORY.create(inverseLhs, generatedTensor.generatedTensor);

//        Transformation[] transformations1 = new Transformation[transformations.length + 2];
//        transformations1[0] = ContractIndices.INSTANCE;
//        System.arraycopy(transformations, 0, transformations1, 1, transformations.length);
//        transformations1[transformations.length + 1] = CalculateNumbers.INSTANCE;
//        final Transformation expandCollect =
//                new ExpandAndCollectTransformation(
//                FullScalarsSplitCriteria.INSTANCE,
//                Indicator.SYMBOL_INDICATOR,
//                transformations1);
//
//        equation.eval(
//                toInverse.asSubstitution(),
//                inverse.asSubstitution(),
//                expandCollect);
        Tensor temp = equation;
        temp = toInverse.transform(temp);
        temp = inverse.transform(temp);
        transformations = ArraysUtils.addAll(new Transformation[]{ContractIndices.INSTANCE}, transformations);
        temp = Expand.expand(temp, transformations);
        temp = ContractIndices.contract(temp);
        equation = (Expression) temp;

        List<Split> rightSplit = new ArrayList<>();

        if (equation.get(1) instanceof Sum)
            for (Tensor summand : equation.get(1))
                rightSplit.add(Split.split(summand));
        else
            rightSplit.add(Split.split(equation.get(1)));
        linearEquations = new ArrayList<>();
        for (Tensor summand : equation.get(0)) {
            Split current = Split.split(summand);
            boolean one = false;
            for (Split split : rightSplit)
                if (TensorUtils.equals(current.summand, split.factor)) {
                    linearEquations.add(ExpressionFactory.FACTORY.create(current.summand, split.factor));
                    one = true;
                    break;
                }
            if (!one)
                linearEquations.add(ExpressionFactory.FACTORY.create(current.summand, Complex.ZERO));
        }
    }

    public Expression inverse() {
        return inverse;
    }

//    private static class Split {
//
//        final Tensor scalar;
//        final Tensor nonScalar;
//
//        public Split(Tensor scalar, Tensor nonScalar) {
//            this.scalar = scalar;
//            this.nonScalar = nonScalar;
//        }
//    }
//    private static Split split(Tensor tensor) {
//        Tensor scalar = null, nonScalar = null;
//        if (tensor instanceof Product) {
//            ProductContent pc = (ProductContent) tensor.getContent();
//            if (pc.getScalarContents().length == 0) {
//                nonScalar = new Product(pc.getRange(0, pc.size()));
//                scalar = new TensorNumber(pc.getFactor());
//            } else {
//                scalar = new Product();
//                if (!pc.getFactor().isOne())
//                    ((Product) scalar).add(new TensorNumber(pc.getFactor()));
//                for (ProductContent content : pc.getScalarContents())
//                    ((Product) scalar).add(content.getRange(0, content.size()));
//                nonScalar = new Product(pc.getNonScalarContent().
//                        getRange(0, pc.getNonScalarContent().size()));
//            }
//        } else {
//            nonScalar = tensor;
//            scalar = TensorNumber.createONE();
//        }
//        return new Split(scalar.equivalent(), nonScalar.equivalent());
//    }
    public void generateMapleFile(String path) {
        try {
            FileOutputStream output = new FileOutputStream(path + "/equations.maple");
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

            file.println("file:=fopen(\"" + path + "/equations.mapleOut\",WRITE);");
            file.append("for k from 1 to " + variables.length + " do\n");
            file.append("fprintf(file,\"%a=%a;\\n\",lhs(Result[1][k]),rhs(Result[1][k]))\n");
            file.append("od;\n");
            file.append("fclose(file);");
        } catch (Exception e) {
        }
    }
}
