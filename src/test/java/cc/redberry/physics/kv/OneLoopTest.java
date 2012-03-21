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
package cc.redberry.physics.kv;

import cc.redberry.core.context.CC;
import cc.redberry.core.context.ToStringMode;
import cc.redberry.core.tensor.*;
import cc.redberry.physics.util.SqrSubs;
import cc.redberry.transformation.CalculateNumbers;
import cc.redberry.transformation.ExpandBrackets;
import cc.redberry.transformation.Transformation;
import cc.redberry.transformation.Transformer;
import cc.redberry.transformation.collect.CollectFactory;
import cc.redberry.transformation.concurrent.EACScalars;
import cc.redberry.transformation.contractions.IndicesContractionsTransformation;
import java.util.Arrays;
import org.junit.Test;
import static cc.redberry.physics.TAssert.*;
import cc.redberry.transformation.Transformations;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class OneLoopTest {
    @Test
    public void RR5() {
        Expression t = new Expression("RR=L*L*(L-1)*HATK^{\\gamma}*DELTA^{\\alpha\\beta}*HATK^{\\mu\\nu}*n_{\\lambda}*((1/20)*R_{\\alpha\\nu}*R^{\\lambda}_{\\gamma\\beta\\mu}+(1/20)*R_{\\alpha\\gamma}*R^{\\lambda}_{\\mu\\beta\\nu}+(1/10)*R_{\\alpha\\beta}*R^{\\lambda}_{\\mu\\gamma\\nu}+(1/20)*R^{\\sigma}_{\\alpha\\nu\\gamma}*R^{\\lambda}_{\\sigma\\beta\\mu}-(1/60)*R^{\\sigma}_{\\mu\\alpha\\nu}*R^{\\lambda}_{\\beta\\sigma\\gamma}+(1/10)*R^{\\sigma}_{\\alpha\\beta\\gamma}*R^{\\lambda}_{\\mu\\sigma\\nu}-(1/12)*R^{\\sigma}_{\\alpha\\beta\\nu}*R^{\\lambda}_{\\mu\\sigma\\gamma})");
        OneLoop loop = new OneLoop();
        t.eval(
                loop.L.asSubstitution(),
                CalculateNumbers.INSTANCE,
                loop.RICCI.asSubstitution(),
                loop.RIMAN.asSubstitution(),
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                CollectFactory.createCollectEqualTerms(),
                CalculateNumbers.INSTANCE);
        System.out.println(t.toString(ToStringMode.UTF8));
    }

    @Test
    public void testINV() {
        OneLoop loop = new OneLoop();
        Expression Kn = new Expression("Kn^{\\alpha\\beta}_{\\gamma\\delta}=K^{\\mu\\nu}^{\\alpha\\beta}_{\\gamma\\delta}*n_{\\mu}*n_{\\nu}");
        Transformation tr = new SqrSubs((SimpleTensor) CC.parse("n_{\\alpha}"));
        Kn.eval(
                loop.MATRIX_K.asSubstitution(),
                loop.P.asSubstitution(),
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                new Transformer(tr),
                OneLoop.KRONECKER_DIMENSION.asSubstitution(),
                CollectFactory.createCollectEqualTerms(),
                CalculateNumbers.INSTANCE,
                EACScalars.getTransformer(),
                CalculateNumbers.INSTANCE);
        Expression One = new Expression("One^{\\alpha\\beta}_{\\rho\\tau}=Kn^{\\alpha\\beta}_{\\gamma\\delta}*KINV^{\\gamma\\delta}_{\\rho\\tau}");
        One.eval(
                Kn.asSubstitution(),
                loop.MATRIX_K_INV.asSubstitution(),
                loop.P.asSubstitution(),
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                OneLoop.KRONECKER_DIMENSION.asSubstitution(),
                new Transformer(tr),
                CollectFactory.createCollectEqualTerms(),
                CalculateNumbers.INSTANCE,
                EACScalars.getTransformer(),
                CalculateNumbers.INSTANCE);
        System.out.println(One);

        Tensor t = CC.parse("1/2*n^{\\nu }*n_{\\nu }*d_{\\rho }^{\\alpha }*d_{\\tau }^{\\beta }");
        System.out.println(tr.transform(t));
    }

    @Test
    public void testInverse() {
        Expression P =
                new Expression("P^{\\alpha\\beta}_{\\mu\\nu} = (1/2)*(d^{\\alpha}_{\\mu}*d^{\\beta}_{\\nu}+d^{\\alpha}_{\\nu}*d^{\\beta}_{\\mu})-"
                + "(1/4)*g_{\\mu\\nu}*g^{\\alpha\\beta}");
        Expression KRONECKER_DIMENSION =
                new Expression("d^{\\alpha}_{\\alpha} = 4");
        Expression MATRIX_K =
                new Expression("K^{\\mu\\nu}^{\\alpha\\beta}_{\\gamma\\delta} = g^{\\mu\\nu}*P^{\\alpha\\beta}_{\\gamma\\delta}+"
                + "(1+2*beta)*((1/4)*(d^{\\mu}_{\\gamma}*g^{\\alpha \\nu}*d^{\\beta}_{\\delta} + d^{\\mu}_{\\delta}*g^{\\alpha \\nu}*d^{\\beta}_{\\gamma}+d^{\\mu}_{\\gamma}*g^{\\beta \\nu}*d^{\\alpha}_{\\delta}+ d^{\\mu}_{\\delta}*g^{\\beta \\nu}*d^{\\alpha}_{\\gamma})+"
                + "(1/4)*(d^{\\nu}_{\\gamma}*g^{\\alpha \\mu}*d^{\\beta}_{\\delta} + d^{\\nu}_{\\delta}*g^{\\alpha \\mu}*d^{\\beta}_{\\gamma}+d^{\\nu}_{\\gamma}*g^{\\beta \\mu}*d^{\\alpha}_{\\delta}+ d^{\\nu}_{\\delta}*g^{\\beta \\mu}*d^{\\alpha}_{\\gamma}) -"
                + "(1/4)*(g_{\\gamma\\delta}*g^{\\mu \\alpha}*g^{\\nu \\beta}+g_{\\gamma\\delta}*g^{\\mu \\beta}*g^{\\nu \\alpha})-"
                + "(1/4)*(g^{\\alpha\\beta}*d^{\\mu}_{\\gamma}*d^{\\nu}_{\\delta}+g^{\\alpha\\beta}*d^{\\mu}_{\\delta}*d^{\\nu}_{\\gamma})+(1/8)*g^{\\mu\\nu}*g_{\\gamma\\delta}*g^{\\alpha\\beta})");

        Transformation tr = new SqrSubs((SimpleTensor) CC.parse("n_{\\alpha}"));
        //Kn
        Expression Kn = new Expression("Kn^{\\alpha\\beta}_{\\gamma\\delta}=K^{\\mu\\nu}^{\\alpha\\beta}_{\\gamma\\delta}*n_{\\mu}*n_{\\nu}");
        Kn.eval(
                MATRIX_K.asSubstitution(),
                P.asSubstitution(),
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                new Transformer(tr),
                KRONECKER_DIMENSION.asSubstitution(),
                CollectFactory.createCollectEqualTerms(),
                CalculateNumbers.INSTANCE,
                EACScalars.getTransformer(),
                CalculateNumbers.INSTANCE);

        //Kn (KV)
        Expression KnKV =
                new Expression("KnKv^{\\alpha\\beta}_{\\gamma\\delta} = P^{\\alpha\\beta}_{\\mu\\nu}*(1/2*d^{\\mu}_{\\rho}*d^{\\nu}_{\\tau}+1/2*d^{\\nu}_{\\rho}*d^{\\mu}_{\\tau}+1/2*(1+2*beta)*(n^{\\mu}*n_{\\rho}*d^{\\nu}_{\\tau}+n^{\\nu}*n_{\\rho}*d^{\\mu}_{\\tau}+n^{\\mu}*n_{\\tau}*d^{\\nu}_{\\rho}+n^{\\nu}*n_{\\tau}*d^{\\mu}_{\\rho}))*P^{\\rho\\tau}_{\\gamma\\delta}");
        KnKV.eval(
                P.asSubstitution(),
                MATRIX_K.asSubstitution(),
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                new Transformer(tr),
                KRONECKER_DIMENSION.asSubstitution(),
                CollectFactory.createCollectEqualTerms(),
                CalculateNumbers.INSTANCE,
                EACScalars.getTransformer(),
                CalculateNumbers.INSTANCE);
        System.out.println(Kn);
        System.out.println(KnKV);

        Tensor diff = new Sum(Kn.right().clone(), new Product(TensorNumber.createMINUSONE(), KnKV.right().clone()));
        diff = new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS).transform(diff);
        diff = CollectFactory.createCollectEqualTerms().transform(diff);
        diff = CalculateNumbers.INSTANCE.transform(diff);
        diff = EACScalars.getTransformer().transform(diff);
        diff = CalculateNumbers.INSTANCE.transform(diff);
        System.out.println("KV - my = " + diff);

        Expression MATRIX_K_INV =
                new Expression("KINV^{\\alpha\\beta}_{\\mu\\nu} = P^{\\alpha\\beta}_{\\mu\\nu}+a*g_{\\mu\\nu}*g^{\\alpha\\beta}+"
                + "(1/4)*b*(n_{\\mu}*n^{\\alpha}*d^{\\beta}_{\\nu}+n_{\\mu}*n^{\\beta}*d^{\\alpha}_{\\nu}+n_{\\nu}*n^{\\alpha}*d^{\\beta}_{\\mu}+n_{\\nu}*n^{\\beta}*d^{\\alpha}_{\\mu})+"
                + "c*(n_{\\mu}*n_{\\nu}*g^{\\alpha\\beta}+n^{\\alpha}*n^{\\beta}*g_{\\mu\\nu})+d*n_{\\mu}*n_{\\nu}*n^{\\alpha}*n^{\\beta}");

        Expression One = new Expression("One^{\\alpha\\beta}_{\\rho\\tau}=KnKv^{\\alpha\\beta}_{\\gamma\\delta}*KINV^{\\gamma\\delta}_{\\rho\\tau}");
        One.eval(
                KnKV.asSubstitution(),
                MATRIX_K_INV.asSubstitution(),
                P.asSubstitution(),
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                KRONECKER_DIMENSION.asSubstitution(),
                new Transformer(tr),
                CollectFactory.createCollectEqualTerms(),
                CalculateNumbers.INSTANCE,
                EACScalars.getTransformer(),
                CalculateNumbers.INSTANCE);
        for (Tensor t : One.right())
            System.out.println("Factor" + Arrays.toString(((ProductContent) t.getContent()).getScalarContents())); //            String s = Arrays.toString(((ProductContent) t.getContent()).getScalarContents());
        //            System.out.println(s.substring(1, s.length() - 1));
//        Expression a = new Expression("a=-(1+2*beta)/(4*(5+6*beta))");
//        Expression b = new Expression("b=-(1+2*beta)/(1+beta)");
//        Expression c = new Expression("c=(1+2*beta)/(5+6*beta)");
//        Expression d = new Expression("d=(1+2*beta)*(1+2*beta)/((1+beta)*(5+6*beta))");
        System.out.println(One);
    }

    @Test
    public void testInverseZeroBeta() {
        Expression P =
                new Expression("P^{\\alpha\\beta}_{\\mu\\nu} = (1/2)*(d^{\\alpha}_{\\mu}*d^{\\beta}_{\\nu}+d^{\\alpha}_{\\nu}*d^{\\beta}_{\\mu})-"
                + "(1/4)*g_{\\mu\\nu}*g^{\\alpha\\beta}");
        Expression KRONECKER_DIMENSION =
                new Expression("d^{\\alpha}_{\\alpha} = 4");
        Expression MATRIX_K =
                new Expression("K^{\\mu\\nu}^{\\alpha\\beta}_{\\gamma\\delta} = g^{\\mu\\nu}*P^{\\alpha\\beta}_{\\gamma\\delta}+"
                + "(1/4)*(d^{\\mu}_{\\gamma}*g^{\\alpha \\nu}*d^{\\beta}_{\\delta} + d^{\\mu}_{\\delta}*g^{\\alpha \\nu}*d^{\\beta}_{\\gamma}+d^{\\mu}_{\\gamma}*g^{\\beta \\nu}*d^{\\alpha}_{\\delta}+ d^{\\mu}_{\\delta}*g^{\\beta \\nu}*d^{\\alpha}_{\\gamma})+"
                + "(1/4)*(d^{\\nu}_{\\gamma}*g^{\\alpha \\mu}*d^{\\beta}_{\\delta} + d^{\\nu}_{\\delta}*g^{\\alpha \\mu}*d^{\\beta}_{\\gamma}+d^{\\nu}_{\\gamma}*g^{\\beta \\mu}*d^{\\alpha}_{\\delta}+ d^{\\nu}_{\\delta}*g^{\\beta \\mu}*d^{\\alpha}_{\\gamma}) -"
                + "(1/4)*(g_{\\gamma\\delta}*g^{\\mu \\alpha}*g^{\\nu \\beta}+g_{\\gamma\\delta}*g^{\\mu \\beta}*g^{\\nu \\alpha})-"
                + "(1/4)*(g^{\\alpha\\beta}*d^{\\mu}_{\\gamma}*d^{\\nu}_{\\delta}+g^{\\alpha\\beta}*d^{\\mu}_{\\delta}*d^{\\nu}_{\\gamma})+(1/8)*g^{\\mu\\nu}*g_{\\gamma\\delta}*g^{\\alpha\\beta}");

        Transformation tr = new SqrSubs((SimpleTensor) CC.parse("n_{\\alpha}"));
        //Kn
        Expression Kn = new Expression("Kn^{\\alpha\\beta}_{\\gamma\\delta}=K^{\\mu\\nu}^{\\alpha\\beta}_{\\gamma\\delta}*n_{\\mu}*n_{\\nu}");
        Kn.eval(
                MATRIX_K.asSubstitution(),
                P.asSubstitution(),
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                new Transformer(tr),
                KRONECKER_DIMENSION.asSubstitution(),
                CollectFactory.createCollectEqualTerms(),
                CalculateNumbers.INSTANCE,
                EACScalars.getTransformer(),
                CalculateNumbers.INSTANCE);

        Expression MATRIX_K_INV =
                new Expression("KINV^{\\alpha\\beta}_{\\mu\\nu} = P^{\\alpha\\beta}_{\\mu\\nu}-(1/20)*g_{\\mu\\nu}*g^{\\alpha\\beta}-"
                + "(1/4)*(n_{\\mu}*n^{\\alpha}*d^{\\beta}_{\\nu}+n_{\\mu}*n^{\\beta}*d^{\\alpha}_{\\nu}+n_{\\nu}*n^{\\alpha}*d^{\\beta}_{\\mu}+n_{\\nu}*n^{\\beta}*d^{\\alpha}_{\\mu})+"
                + "(1/5)*(n_{\\mu}*n_{\\nu}*g^{\\alpha\\beta}+n^{\\alpha}*n^{\\beta}*g_{\\mu\\nu})+(1/5)*n_{\\mu}*n_{\\nu}*n^{\\alpha}*n^{\\beta}");

        Expression One = new Expression("One^{\\alpha\\beta}_{\\rho\\tau}=Kn^{\\alpha\\beta}_{\\gamma\\delta}*KINV^{\\gamma\\delta}_{\\rho\\tau}");
        One.eval(
                Kn.asSubstitution(),
                MATRIX_K_INV.asSubstitution(),
                P.asSubstitution(),
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                KRONECKER_DIMENSION.asSubstitution(),
                new Transformer(tr),
                CollectFactory.createCollectEqualTerms(),
                CalculateNumbers.INSTANCE,
                EACScalars.getTransformer(),
                CalculateNumbers.INSTANCE);
        Tensor expected = CC.parse("P^{\\alpha\\beta}_{\\rho\\tau}");
        expected = P.asSubstitution().transform(expected);
        expected = Transformations.expandBrackets(expected);
        expected = Transformations.calculateNumbers(expected);
        assertParity(expected, One.right());
    }

    @Test
    public void testInverseZeroBeta1() {
        Expression P =
                new Expression("P^{\\alpha\\beta}_{\\mu\\nu} = (1/2)*(d^{\\alpha}_{\\mu}*d^{\\beta}_{\\nu}+d^{\\alpha}_{\\nu}*d^{\\beta}_{\\mu})-"
                + "(1/4)*g_{\\mu\\nu}*g^{\\alpha\\beta}");
        Expression KRONECKER_DIMENSION =
                new Expression("d^{\\alpha}_{\\alpha} = 4");
        Expression MATRIX_K =
                new Expression("K^{\\mu\\nu}^{\\alpha\\beta}_{\\gamma\\delta} = g^{\\mu\\nu}*P^{\\alpha\\beta}_{\\gamma\\delta}+"
                + "(1/4)*(d^{\\mu}_{\\gamma}*g^{\\alpha \\nu}*d^{\\beta}_{\\delta} + d^{\\mu}_{\\delta}*g^{\\alpha \\nu}*d^{\\beta}_{\\gamma}+d^{\\mu}_{\\gamma}*g^{\\beta \\nu}*d^{\\alpha}_{\\delta}+ d^{\\mu}_{\\delta}*g^{\\beta \\nu}*d^{\\alpha}_{\\gamma})+"
                + "(1/4)*(d^{\\nu}_{\\gamma}*g^{\\alpha \\mu}*d^{\\beta}_{\\delta} + d^{\\nu}_{\\delta}*g^{\\alpha \\mu}*d^{\\beta}_{\\gamma}+d^{\\nu}_{\\gamma}*g^{\\beta \\mu}*d^{\\alpha}_{\\delta}+ d^{\\nu}_{\\delta}*g^{\\beta \\mu}*d^{\\alpha}_{\\gamma}) -"
                + "(1/4)*(g_{\\gamma\\delta}*g^{\\mu \\alpha}*g^{\\nu \\beta}+g_{\\gamma\\delta}*g^{\\mu \\beta}*g^{\\nu \\alpha})-"
                + "(1/4)*(g^{\\alpha\\beta}*d^{\\mu}_{\\gamma}*d^{\\nu}_{\\delta}+g^{\\alpha\\beta}*d^{\\mu}_{\\delta}*d^{\\nu}_{\\gamma})+(1/8)*g^{\\mu\\nu}*g_{\\gamma\\delta}*g^{\\alpha\\beta}");

        Transformation tr = new SqrSubs((SimpleTensor) CC.parse("n_{\\alpha}"));
        Expression MATRIX_K_INV =
                new Expression("KINV^{\\alpha\\beta}_{\\mu\\nu} = P^{\\alpha\\beta}_{\\mu\\nu}-(1/20)*g_{\\mu\\nu}*g^{\\alpha\\beta}-"
                + "(1/4)*(n_{\\mu}*n^{\\alpha}*d^{\\beta}_{\\nu}+n_{\\mu}*n^{\\beta}*d^{\\alpha}_{\\nu}+n_{\\nu}*n^{\\alpha}*d^{\\beta}_{\\mu}+n_{\\nu}*n^{\\beta}*d^{\\alpha}_{\\mu})+"
                + "(1/5)*(n_{\\mu}*n_{\\nu}*g^{\\alpha\\beta}+n^{\\alpha}*n^{\\beta}*g_{\\mu\\nu})+(1/5)*n_{\\mu}*n_{\\nu}*n^{\\alpha}*n^{\\beta}");

        Expression One = new Expression("One^{\\mu\\nu}^{\\alpha\\beta}_{\\rho\\tau}=K^{\\mu\\nu}^{\\alpha\\beta}_{\\gamma\\delta}*KINV^{\\gamma\\delta}_{\\rho\\tau}");
        One.eval(
                MATRIX_K.asSubstitution(),
                MATRIX_K_INV.asSubstitution(),
                P.asSubstitution(),
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                KRONECKER_DIMENSION.asSubstitution(),
                new Transformer(tr),
                CollectFactory.createCollectEqualTerms(),
                CalculateNumbers.INSTANCE,
                EACScalars.getTransformer(),
                CalculateNumbers.INSTANCE);
        Expression OneN = new Expression("One^{\\alpha\\beta}_{\\rho\\tau} = One^{\\mu\\nu}^{\\alpha\\beta}_{\\rho\\tau}*n_{\\mu}*n_{\\nu}");
        OneN.eval(
                One.asSubstitution(),
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                KRONECKER_DIMENSION.asSubstitution(),
                new Transformer(tr),
                CollectFactory.createCollectEqualTerms(),
                CalculateNumbers.INSTANCE,
                EACScalars.getTransformer(),
                CalculateNumbers.INSTANCE);

        Tensor expected = CC.parse("P^{\\alpha\\beta}_{\\rho\\tau}");
        expected = P.asSubstitution().transform(expected);
        expected = Transformations.expandBrackets(expected);
        expected = Transformations.calculateNumbers(expected);
        assertParity(expected, OneN.right());
    }
}
