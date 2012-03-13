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
package cc.redberry.physics.kv;

import cc.redberry.core.context.CC;
import cc.redberry.core.context.ToStringMode;
import cc.redberry.core.tensor.*;
import cc.redberry.core.tensor.iterators.TensorFirstTreeIterator;
import cc.redberry.core.tensor.iterators.TensorTreeIterator;
import cc.redberry.core.tensor.testing.TTest;
import cc.redberry.core.utils.Indicator;
import cc.redberry.transformation.*;
import cc.redberry.transformation.collect.CollectFactory;
import cc.redberry.transformation.collect.CollectInputPortImpl;
import cc.redberry.transformation.collect.EqualsSplitCriteria;
import cc.redberry.transformation.collect.ScalarsSplitCriteria;
import cc.redberry.transformation.concurrent.ExpandAndCollectTransformation;
import cc.redberry.transformation.contractions.IndicesContractionsTransformation;
import java.util.ArrayList;
import java.util.List;
import org.junit.Ignore;
import org.junit.Test;
import cc.redberry.physics.util.IndicesFactoryUtil;
import static cc.redberry.physics.TAssert.*;
import static cc.redberry.physics.util.IndicesFactoryUtil.*;
import cc.redberry.physics.util.SqrSubs;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 * @author Konstantin Kiselev
 */
public class OneLoop1Test {
    public OneLoop1Test() {
    }

    @Test
    public void evalRR() {
        OneLoop loop1 = new OneLoop(OneLoop.EVAL.INITIALIZE);
        loop1.evalRR();
        System.out.println(loop1.RR.toString(ToStringMode.UTF8));
    }

    @Test
    public void evalHATK() {
        OneLoop loop1 = new OneLoop(OneLoop.EVAL.INITIALIZE);
        loop1.evalHatK();
        System.out.println(loop1.HATK_1.toString(ToStringMode.UTF8));
    }

    @Test
    public void evalDeltas() {
        OneLoop loop1 = new OneLoop(OneLoop.EVAL.INITIALIZE);
        loop1.evalDeltas();
        Sum s = (Sum) loop1.DELTA_3.right();
        loop1.evalHatK();
        CollectInputPortImpl cip = new CollectInputPortImpl(EqualsSplitCriteria.INSTANCE);
        System.out.println(s.size());
        Product p = (Product) s.getElements().get(6);
        for (Expression hatk : loop1.HATKs)
            hatk.asSubstitution().transform(p);

        Tensor p1 = new Product(p.getElements().get(1).clone(), p.getElements().get(2).clone());
        p1 = ExpandAndCollectTransformation.EXPAND_AND_COLLECT.transform(p1);
        Tensor p2 = new Product(p1, p.getElements().get(3).clone());
        p2 = ExpandAndCollectTransformation.EXPAND_AND_COLLECT.transform(p2);
        System.out.println(p2.toString(ToStringMode.UTF8));
//        for (int i = 0; i < s.size(); ++i) {
//            Tensor t = s.getElements().get(i);
//
//            for (Expression hatk : loop1.HATKs)
//                t = hatk.asSubstitution().transform(t);
//            OutputPortUnsafe<Tensor> bracketsOutput = ExpandBracketsOutput.create(t, Indicator.SYMBOL_INDICATOR);
//            Tensor cur;
//            while ((cur = bracketsOutput.take()) != null) {
//                cur = Transformations.contractMetrics(cur);
//                cur = loop1.KRONECKER_DIMENSION.asSubstitution().transform(cur);
//                cur = Transformations.calculateNumbers(cur);
//                cip.put(cur);
//            }
//            System.out.println(i);
//        }
    }

    @Ignore
    @Test
    public void test() {
        OneLoop loop1 = new OneLoop(OneLoop.EVAL.INITIALIZE);
        loop1.evalRR();
        assertIndices(loop1.RR);
        loop1.evalDeltas();
        System.out.println(loop1.DELTA_3.toString(ToStringMode.UTF8));
        System.out.println(loop1.DELTA_2.toString(ToStringMode.UTF8));
        System.out.println(loop1.DELTA_1.toString(ToStringMode.UTF8));
        System.out.println(loop1.DELTA_4.toString(ToStringMode.UTF8));
        assertIndices(loop1.DELTAs);
//        loop1.DELTA_3.asSubstitution();
        loop1.DELTA_4.asSubstitution();
    }

    @Ignore
    @Test
    public void testHATKs() {
        OneLoop loop1 = new OneLoop(OneLoop.EVAL.EVAL_HATK);
        System.out.println(loop1.HATK_1.toString(ToStringMode.UTF8));
        System.out.println(loop1.HATK_2.toString(ToStringMode.UTF8));
        System.out.println(loop1.DELTA_1.toString(ToStringMode.UTF8));
        System.out.println(loop1.DELTA_2.toString(ToStringMode.UTF8));
    }

    @Ignore
    @Test
    public void testAll() {
        OneLoop loop1 = new OneLoop(OneLoop.EVAL.EVAL_ALL);
        System.out.println(loop1.RR.toString(ToStringMode.UTF8));
        System.out.println(((MultiTensor) loop1.RR.right()).size());
        assertIndices(loop1.RR);
    }

    @Ignore
    @Test
    public void testDelta2() {
        OneLoop loop1 = new OneLoop(OneLoop.EVAL.EVAL_HATK);
        IndicesInsertion indicesInsertion = new IndicesInsertion(loop1.matricesIndicator, IndicesFactoryUtil.createIndices(loop1.DELTAs, "^{\\mu\\nu}_{\\alpha\\beta}"));

        loop1.DELTA_2.eval(
                indicesInsertion,
                loop1.L.asSubstitution(),
                CalculateNumbers.INSTANCE);

        System.out.println("HATK subs ... ");
        loop1.DELTA_2.eval(
                loop1.HATK_1.asSubstitution(),
                loop1.HATK_2.asSubstitution(),
                loop1.HATK_3.asSubstitution(),
                loop1.HATK_4.asSubstitution());
        System.out.println("HATK subs ... done");

        System.out.println("Expand brackets ...");
        loop1.DELTA_2.eval(
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS));
        System.out.println("Expand brackets ... done");

        System.out.println("Indices contractions ...");
        loop1.DELTA_2.eval(
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                loop1.KRONECKER_DIMENSION.asSubstitution(),
                CalculateNumbers.INSTANCE);
        System.out.println("Indices contractions ... done");
        assertIndices(loop1.DELTA_2.right());
        System.out.println("Collecting " + ((MultiTensor) loop1.DELTA_2.right()).size() + " elements ...");

        loop1.DELTA_2.eval(
                CollectFactory.createCollectEqualTerms());
        System.out.println("Collecting ... done. The new size: " + ((MultiTensor) loop1.DELTA_2.right()).size());

        System.out.println("Scalars ...");
        loop1.DELTA_2.eval(
                CalculateNumbers.INSTANCE,
                CollectFactory.createCollectAllScalars(),
                CalculateNumbers.INSTANCE);
        System.out.println("Scalars ... done");

        System.out.println(loop1.DELTA_2.toString(ToStringMode.UTF8));
        Tensor rhs = loop1.DELTA_2.right();
        TensorTreeIterator iterator = new TensorFirstTreeIterator(rhs);
        while (iterator.hasNext())
            if (TTest.testIsSymbol(iterator.next()))
                iterator.remove();
        System.out.println(rhs.toString(ToStringMode.UTF8));
    }

    @Ignore
    @Test
    public void test1() {
        OneLoop loop1 = new OneLoop(OneLoop.EVAL.INITIALIZE);
        loop1.evalRR();
        loop1.evalDeltas();
        for (Expression delta : loop1.DELTAs)
            loop1.RR.eval(delta.asSubstitution());
        loop1.RR.eval(
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                loop1.KRONECKER_DIMENSION.asSubstitution(),
                CalculateNumbers.INSTANCE,
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                CollectFactory.createCollectEqualTerms(),
                CollectFactory.createCollectScalars(),
                CalculateNumbers.INSTANCE);

        Transformation indicesInsertion;
        indicesInsertion = new IndicesInsertion(loop1.matricesIndicator, createIndices(loop1.HATKs, "^{\\mu\\nu}_{\\alpha\\beta}"));
        for (Expression hatK : loop1.HATKs) {
            hatK.eval(indicesInsertion);
            loop1.RR.eval(hatK.asSubstitution());
        }


//        loop1.MATRIX_K.eval(
//                loop1.P.asSubstitution(),
//                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
//                CollectFactory.createCollectEqualTerms(),
//                CollectFactory.createCollectEqualTerms());
        System.out.println("matrix K evaluating");

        loop1.MATRIX_K.eval(loop1.P.asSubstitution(),
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                CollectFactory.createCollectEqualTerms(),
                CalculateNumbers.INSTANCE);
        loop1.MATRIX_K_INV.eval(loop1.P.asSubstitution(),
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                CollectFactory.createCollectEqualTerms(),
                CalculateNumbers.INSTANCE);
        System.out.println("done");
        loop1.RR.eval(IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC);
        System.out.println(((MultiTensor) loop1.RR.right()).getElements().get(0).toString(ToStringMode.UTF8));
        loop1.RR.eval(loop1.MATRIX_K_INV.asSubstitution());//, CollectFactory.createCollectEqualTerms(), CalculateNumbers.INSTANCE, CollectFactory.createCollectAllScalars());
        loop1.RR.eval(loop1.MATRIX_K.asSubstitution());


        Tensor rhs = ((MultiTensor) loop1.RR.right()).getElements().get(0).clone();
        System.out.println(" Evaluating RR ");
        loop1.RR.eval(
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                loop1.KRONECKER_DIMENSION.asSubstitution());

        TensorTreeIterator treeIterator = new TensorFirstTreeIterator(rhs);
        Tensor current;
        while (treeIterator.hasNext()) {
            current = treeIterator.next();
            if (current instanceof Sum && !TTest.testIsSymbol(current)) {
                current = CollectFactory.createCollectEqualTerms().transform(current.clone());
                current = CalculateNumbers.INSTANCE.transform(current);
                treeIterator.set(current);
            }
        }
        List<Sum> sums = new ArrayList<>();
        TensorIterator iterator = rhs.iterator();
        while (iterator.hasNext()) {
            current = iterator.next();
            if (current instanceof Sum && !TTest.testIsSymbol(current)) {
                sums.add((Sum) current.clone());
                iterator.remove();
            }
        }


        System.out.println(" counting terms: ... ");
        long start = System.currentTimeMillis();
        Transformation ec = new ExpandAndCollectTransformation(
                Indicator.SYMBOL_INDICATOR,
                new Transformation[]{IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    loop1.KRONECKER_DIMENSION.asSubstitution(),
                    CalculateNumbers.INSTANCE});
        for (Sum sum : sums) {
            rhs = new Product(rhs.equivalent(), sum);
            rhs = ec.transform(rhs.clone());
            System.out.println("+");
        }
        rhs = loop1.KRONECKER_DIMENSION.asSubstitution().transform(rhs);
        rhs = Transformations.calculateNumbers(rhs);
        long stop = System.currentTimeMillis();
        System.out.println("... . Evalution time = " + (stop - start) + "ms");
        System.out.println(rhs.toString(ToStringMode.UTF8));
//        System.out.println(fS.transform(rhs).toString(ToStringMode.UTF8));

//        Tensor rhs_1 = ((MultiTensor) rhs).getElements().get(0);
//
//        System.out.println("contracting  ");
//        rhs_1 = IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC.transform(rhs_1);
//        rhs_1 = loop1.KRONECKER_DIMENSION.asSubstitution().transform(rhs_1);
//        rhs_1 = CalculateNumbers.INSTANCE.transform(rhs_1);
//        System.out.println("done");
//        System.out.println("expanding  ");
//        rhs_1 = new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS).transform(rhs_1);
//        System.out.println("done");
//        System.out.println("contracting  ");
//        rhs_1 = IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC.transform(rhs_1);
//        rhs_1 = loop1.KRONECKER_DIMENSION.asSubstitution().transform(rhs_1);
//        rhs_1 = CalculateNumbers.INSTANCE.transform(rhs_1);
//        System.out.println("done");
//        System.out.println(rhs_1.toString(ToStringMode.UTF8));
//
//        int oldSize = ((MultiTensor) rhs_1).size();
//        TensorIterator iterator = rhs_1.iterator();
//        Tensor c;
//        int count;
//        while (iterator.hasNext()) {
//            c = iterator.next();
//            count = 0;
//            for (Tensor n : c)
//                if (TTest.testEqualstensorStructure(n, CC.parse("n_{\\mu}")))
//                    count++;
//            if (count % 2 != 0)
//                iterator.remove();
//        }
//        int newSize = ((MultiTensor) rhs_1).size();
//        System.out.println("Old size: " + oldSize + ", new sizeL " + newSize + ".  Collecting");
//        Tensor nT = ((MultiTensor) rhs_1).getElements().get(0);
////        count = 0;
////        for (Tensor t : rhs_1) {
////            if (count > 10)
////                break;
////            ((Sum) nT).add(t);
////        }
//        System.out.println(nT.toString(ToStringMode.UTF8));
//        nT = loop1.MATRIX_K.asSubstitution().transform(nT);
//        nT = Transformations.expandBrackets(nT);
//        nT = Transformations.contractMetrics(nT);
//        nT = loop1.KRONECKER_DIMENSION.asSubstitution().transform(nT);
//        nT = Transformations.calculateNumbers(nT);
//        System.out.println("Collecting");
//
//        nT = CollectFactory.createCollectEqualTerms().transform(nT);
//        System.out.println(((Sum) nT).size());
//        System.out.println(nT.toString(ToStringMode.UTF8));
    }

//    @Ignore
    @Test
    public void test2() {
        OneLoop loop1 = new OneLoop(OneLoop.EVAL.INITIALIZE);
        loop1.evalRR();
        loop1.evalHatK();
        loop1.evalDeltas();
        for (Expression delta : loop1.DELTAs)
            loop1.RR.eval(delta.asSubstitution());
        for (Expression hatK : loop1.HATKs)
            loop1.RR.eval(hatK.asSubstitution());
        loop1.RR.eval(
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                loop1.KRONECKER_DIMENSION.asSubstitution());

        Tensor rhs = ((MultiTensor) loop1.RR.right()).getElements().get(0).clone();
//        System.out.println(rhs.toString(ToStringMode.UTF8));
        System.out.println(" Evaluating RR ");

        TensorTreeIterator treeIterator = new TensorFirstTreeIterator(rhs);
        Tensor current;
//        while (treeIterator.hasNext()) {
//            current = treeIterator.next();
//            if (current instanceof Sum && !TTest.testIsSymbol(current)) {
//                current = CollectFactory.createCollectEqualTerms().transform(current.clone());
//                current = CalculateNumbers.INSTANCE.transform(current);
//                treeIterator.set(current);
//            }
//        }
        List<Sum> sums = new ArrayList<>();
        TensorIterator iterator = rhs.iterator();
        while (iterator.hasNext()) {
            current = iterator.next();
            if (current instanceof Sum && !TTest.testIsSymbol(current)) {
                sums.add((Sum) current.clone());
                iterator.remove();
            }
        }

        System.out.println(" counting terms: ... Sums number " + sums.size());
        Transformation tr = new SqrSubs((SimpleTensor) CC.parse("n_{\\alpha}"));

        long start = System.currentTimeMillis();
        Transformation ec = new ExpandAndCollectTransformation(
                EqualsSplitCriteria.INSTANCE,
                Indicator.SYMBOL_INDICATOR,
                new Transformation[]{
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    loop1.KRONECKER_DIMENSION.asSubstitution(),
                    tr,
                    CalculateNumbers.INSTANCE});
        Product firstPart = new Product(sums.get(1), sums.get(2));
        System.out.println(sums.get(1).size());
        System.out.println(sums.get(2).size());
        Tensor fst = ec.transform(firstPart);
        System.out.println("+");
//        Tensor scn = ec.transform(secondPart);
        System.out.println(((MultiTensor) fst).size());
//        System.out.println(((MultiTensor) scn).size());
        Product fin = new Product();
        fin.add(rhs);
        fin.add(fst);
        fin.add(sums.get(0));
//        
        Transformation sP = new Transformation() {
            @Override
            public Tensor transform(Tensor tensor) {
                if (!(tensor instanceof Product))
                    throw new RuntimeException();
                Product s = new Product();
                for (Tensor t : tensor)
                    if (TTest.testIsSymbol(t))
                        s.add(t);
                if (s.isEmpty())
                    return TensorNumber.createONE();
                else
                    return s.equivalent();
            }
        };
        Transformation ec1 = new ExpandAndCollectTransformation(
                ScalarsSplitCriteria.INSTANCE,
                Indicator.SYMBOL_INDICATOR,
                new Transformation[]{
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    loop1.KRONECKER_DIMENSION.asSubstitution(),
                    CalculateNumbers.INSTANCE,
                    sP});

        Tensor AAA = ec1.transform(rhs);
        System.out.println(((MultiTensor) AAA).size());
//        System.out.println(AAA.toString(ToStringMode.UTF8));
//        for (Sum sum : sums) {
//            rhs = new Product(rhs.equivalent(), sum);
//            rhs = ec.transform(rhs.clone());
//            System.out.println("+");
//        }

//        Transformation ec2 = new ExpandAndCollectTransformation(
//                ScalarsSplitCriteria.INSTANCE,
//                Indicator.FALSE_INDICATOR,
//                new Transformation[]{
//                    CalculateNumbers.INSTANCE});
//        Tensor Res = ec2.transform(AAA);
//        System.out.println(((MultiTensor) Res).size());
//        System.out.println(Res);
//        long stop = System.currentTimeMillis();
//        System.out.println("... . Evalution time = " + (stop - start) + "ms");
//        System.out.println(rhs.toString(ToStringMode.UTF8));
    }

    @Ignore
    @Test
    public void test2_2() {
        OneLoop loop1 = new OneLoop(OneLoop.EVAL.INITIALIZE);
        loop1.evalRR();
        loop1.evalHatK();
        loop1.evalDeltas();
        for (Expression delta : loop1.DELTAs)
            loop1.RR.eval(delta.asSubstitution());
        System.out.println(loop1.RR.toString(ToStringMode.UTF8));
        for (Expression hatK : loop1.HATKs)
            loop1.RR.eval(hatK.asSubstitution());
        loop1.RR.eval(
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                loop1.KRONECKER_DIMENSION.asSubstitution());

        int sC = 0, tC = 0;;
        for (Tensor t : loop1.RR.right()) {
            sC = 0;
            for (Tensor m : t)
                if (m instanceof Sum && !TTest.testIsSymbol(m))
                    sC++;
//            if (sC != 4)
            System.out.println(sC);
            tC++;
        }
        System.out.println(tC);
    }
    
    @Test
    public void testAsasda(){
        Tensor t = CC.parse("1319303/30720*LAMBDA*LAMBDA+20279/240*LAMBDA*LAMBDA*beta+1710821/3840*LAMBDA*LAMBDA*beta*beta+249437/5120*LAMBDA*LAMBDA*b+899449/5120*LAMBDA*LAMBDA*b*beta+298541/320*LAMBDA*LAMBDA*b*beta*beta+1963047/81920*LAMBDA*LAMBDA*b*b+72781753/46080*LAMBDA*LAMBDA*b*b*beta+6613713/10240*LAMBDA*LAMBDA*b*b*beta*beta+215963873/576*LAMBDA*LAMBDA*beta*beta*beta+2162157/640*LAMBDA*LAMBDA*b*beta*beta*beta+2621613/640*LAMBDA*LAMBDA*beta*beta*beta*beta+150879329/2880*LAMBDA*LAMBDA*b*beta*beta*beta*beta+5248947/2560*LAMBDA*LAMBDA*b*beta*beta*b*beta+(-77456989/9216)*LAMBDA*LAMBDA*b*beta*beta*b*beta*beta+282417/80*LAMBDA*LAMBDA*beta*beta*beta*beta*beta+1481397/320*LAMBDA*LAMBDA*b*beta*beta*beta*beta*beta+1583397/640*LAMBDA*LAMBDA*b*beta*beta*b*beta*beta*beta+896877/163840*LAMBDA*LAMBDA*b*b*b+282600331/368640*LAMBDA*LAMBDA*b*b*beta*b+18543097/2880*LAMBDA*LAMBDA*b*beta*beta*b*beta*b+1902681/2048*LAMBDA*LAMBDA*b*beta*beta*b*beta*beta*b+1770417/2560*LAMBDA*LAMBDA*b*beta*beta*b*beta*beta*b*beta+(-2097152/45)*LAMBDA*LAMBDA*b*beta*beta*c+0*LAMBDA*LAMBDA*b*b*beta*beta*c+2097152/45*LAMBDA*LAMBDA*b*beta*beta*b*beta*beta*c+(-535889171/23040)*LAMBDA*LAMBDA*d+159593/960*LAMBDA*LAMBDA*beta*d+119227/192*LAMBDA*LAMBDA*beta*beta*d+385913/10240*LAMBDA*LAMBDA*b*d+69077353/11520*LAMBDA*LAMBDA*b*beta*d+1028787/1280*LAMBDA*LAMBDA*b*beta*beta*d+454023/40960*LAMBDA*LAMBDA*b*b*d+377943/5120*LAMBDA*LAMBDA*b*b*beta*d+(-248042383/46080)*LAMBDA*LAMBDA*b*b*beta*beta*d+27116251/144*LAMBDA*LAMBDA*beta*beta*beta*d+823331/320*LAMBDA*LAMBDA*b*beta*beta*beta*d+491613/160*LAMBDA*LAMBDA*beta*beta*beta*beta*d+2668473/640*LAMBDA*LAMBDA*b*beta*beta*beta*beta*d+284891821/11520*LAMBDA*LAMBDA*b*beta*beta*b*beta*d+(-1021328209/23040)*LAMBDA*LAMBDA*b*beta*beta*b*beta*beta*d+2397*LAMBDA*LAMBDA*beta*beta*beta*beta*d*beta+129453/40*LAMBDA*LAMBDA*b*beta*beta*beta*beta*d*beta+279507/160*LAMBDA*LAMBDA*b*beta*beta*b*beta*beta*d*beta+306510883/184320*LAMBDA*LAMBDA*b*b*b*beta*beta+262144/45*LAMBDA*LAMBDA*c*beta*b*beta*b*b+(-1048576/45)*LAMBDA*LAMBDA*c*beta*b*beta*d*b+87303/5120*LAMBDA*LAMBDA*d*d+16707/320*LAMBDA*LAMBDA*beta*d*d+136131569/5760*LAMBDA*LAMBDA*beta*d*beta*d+54068423/288*LAMBDA*LAMBDA*beta*d*beta*d*beta+(-1072907281/92160)*LAMBDA*LAMBDA*b*d*d+269713069/23040*LAMBDA*LAMBDA*beta*b*d*d+273580171/11520*LAMBDA*LAMBDA*beta*b*beta*d*d+68269999/720*LAMBDA*LAMBDA*beta*b*beta*d*beta*d+736407/320*LAMBDA*LAMBDA*beta*beta*d*beta*d*beta+552104663/5760*LAMBDA*LAMBDA*beta*b*beta*d*beta*d*beta+1761/512*LAMBDA*LAMBDA*d*d*d+5037/320*LAMBDA*LAMBDA*beta*d*d*d+46617/320*LAMBDA*LAMBDA*beta*d*beta*d*d+44811/80*LAMBDA*LAMBDA*beta*d*beta*d*beta*d+148653/160*LAMBDA*LAMBDA*beta*d*beta*d*beta*d*beta+549741893/23592960*LAMBDA*LAMBDA*b*b*b*b+473499/81920*LAMBDA*LAMBDA*beta*b*b*b*b+(-438880433/2949120)*LAMBDA*LAMBDA*beta*b*beta*b*b*b+7885701/81920*LAMBDA*LAMBDA*beta*b*beta*b*beta*b*b+538932443/1474560*LAMBDA*LAMBDA*d*b*b*b+137852603/184320*LAMBDA*LAMBDA*beta*d*b*b*b+2702601/20480*LAMBDA*LAMBDA*beta*d*beta*b*b*b+(-518161127/46080)*LAMBDA*LAMBDA*beta*d*beta*b*beta*b*b+(-324763259/1474560)*LAMBDA*LAMBDA*beta*b*beta*b*beta*b*beta*b+6381027/10240*LAMBDA*LAMBDA*beta*d*beta*b*beta*b*beta*b+(-248655499/184320)*LAMBDA*LAMBDA*beta*b*beta*b*beta*b*beta*b*beta+300843/640*LAMBDA*LAMBDA*beta*d*beta*b*beta*b*beta*b*beta+72441/40*LAMBDA*LAMBDA*beta*d*beta*beta*d*beta*beta+271353319/1440*LAMBDA*LAMBDA*beta*d*beta*b*beta*d*beta*beta+538118231/737280*LAMBDA*LAMBDA*d*b*d*b+(-13278095/4608)*LAMBDA*LAMBDA*beta*d*b*d*b+57996437/18432*LAMBDA*LAMBDA*beta*d*beta*b*d*b+1892079/2560*LAMBDA*LAMBDA*beta*d*beta*b*beta*d*b+588990119/46080*LAMBDA*LAMBDA*beta*d*beta*b*beta*d*beta*b+(-532035617/5760)*LAMBDA*LAMBDA*beta*d*beta*b*beta*d*beta*b*beta+2079/2048*LAMBDA*LAMBDA*d*b*d*d+23823/1280*LAMBDA*LAMBDA*beta*d*b*d*d+194589/1280*LAMBDA*LAMBDA*beta*d*beta*b*d*d+538294649/2880*LAMBDA*LAMBDA*beta*d*beta*b*beta*d*d+541234787/5760*LAMBDA*LAMBDA*beta*d*beta*b*beta*d*beta*d+7047/10*LAMBDA*LAMBDA*beta*d*beta*beta*d*beta*d*beta+22221/40*LAMBDA*LAMBDA*beta*d*beta*b*beta*d*beta*d*beta+585/2048*LAMBDA*LAMBDA*d*d*d*d+117/64*LAMBDA*LAMBDA*beta*d*d*d*d+24399/1280*LAMBDA*LAMBDA*beta*d*beta*d*d*d+22869/320*LAMBDA*LAMBDA*beta*d*beta*d*beta*d*d+74061/640*LAMBDA*LAMBDA*beta*d*beta*d*beta*d*beta*d+1377/16*LAMBDA*LAMBDA*beta*d*beta*d*beta*d*beta*d*beta+4096/15*LAMBDA*LAMBDA*c*beta*b*b*b*b+(-106496/45)*LAMBDA*LAMBDA*c*beta*b*beta*b*b*b+0*LAMBDA*LAMBDA*c*beta*b*beta*b*beta*b*beta*b+5706/5*LAMBDA*LAMBDA*beta*beta*beta*beta*beta*beta+11349/8*LAMBDA*LAMBDA*beta*beta*b*beta*beta*beta*beta+117243/160*LAMBDA*LAMBDA*beta*beta*b*beta*beta*b*beta*beta+64791/320*LAMBDA*LAMBDA*beta*beta*b*beta*b*beta*b*beta*beta+7227/10*LAMBDA*LAMBDA*beta*beta*beta*d*beta*beta*beta+9711/10*LAMBDA*LAMBDA*beta*beta*b*beta*d*beta*beta*beta+67484461/720*LAMBDA*LAMBDA*beta*beta*b*beta*d*beta*b*beta*beta+40743/1280*LAMBDA*LAMBDA*beta*b*beta*b*beta*b*beta*b*beta*beta+11169/80*LAMBDA*LAMBDA*beta*d*beta*b*beta*b*beta*b*beta*beta+5373/10*LAMBDA*LAMBDA*beta*d*beta*beta*d*beta*beta*beta+23787/40*LAMBDA*LAMBDA*beta*d*beta*b*beta*d*beta*beta*beta+38961/160*LAMBDA*LAMBDA*beta*d*beta*b*beta*d*beta*b*beta*beta+405/2*LAMBDA*LAMBDA*beta*d*beta*beta*d*beta*d*beta*beta+3159/20*LAMBDA*LAMBDA*beta*d*beta*b*beta*d*beta*d*beta*beta+243/10*LAMBDA*LAMBDA*beta*d*beta*d*beta*d*beta*d*beta*beta+54999/2621440*LAMBDA*LAMBDA*b*b*b*b*b+567081/1310720*LAMBDA*LAMBDA*beta*b*b*b*b*b+142110611/2949120*LAMBDA*LAMBDA*beta*b*beta*b*b*b*b+279824623/1474560*LAMBDA*LAMBDA*beta*b*beta*b*beta*b*b*b+(-66965251/737280)*LAMBDA*LAMBDA*d*b*b*b*b+29781/8192*LAMBDA*LAMBDA*beta*d*b*b*b*b+175095/8192*LAMBDA*LAMBDA*beta*d*beta*b*b*b*b+151551/2560*LAMBDA*LAMBDA*beta*d*beta*b*beta*b*b*b+274150573/491520*LAMBDA*LAMBDA*beta*b*beta*b*beta*b*beta*b*b+(-15789907/11520)*LAMBDA*LAMBDA*beta*d*beta*b*beta*b*beta*b*b+145557/16384*LAMBDA*LAMBDA*beta*b*beta*b*beta*b*beta*b*beta*b+35004899/23040*LAMBDA*LAMBDA*beta*d*beta*b*beta*b*beta*b*beta*b+42183/81920*LAMBDA*LAMBDA*d*b*d*b*b+271987873/368640*LAMBDA*LAMBDA*beta*d*b*d*b*b+286083/5120*LAMBDA*LAMBDA*beta*d*beta*b*d*b*b+55084811/9216*LAMBDA*LAMBDA*beta*d*beta*b*beta*d*b*b+1097091/5120*LAMBDA*LAMBDA*beta*d*beta*b*beta*d*beta*b*b+271963249/23040*LAMBDA*LAMBDA*beta*d*beta*b*beta*d*beta*b*beta*b+(-67098739/23040)*LAMBDA*LAMBDA*d*b*d*d*b+1179/128*LAMBDA*LAMBDA*beta*d*b*d*d*b+139941/2560*LAMBDA*LAMBDA*beta*d*beta*b*d*d*b+47511/320*LAMBDA*LAMBDA*beta*d*beta*b*beta*d*d*b+(-66510517/2880)*LAMBDA*LAMBDA*beta*d*beta*b*beta*d*beta*d*b+11691/80*LAMBDA*LAMBDA*beta*d*beta*b*beta*d*beta*d*beta*b+135/2048*LAMBDA*LAMBDA*d*d*d*d*b+2457/1024*LAMBDA*LAMBDA*beta*d*d*d*d*b+20169/1280*LAMBDA*LAMBDA*beta*d*beta*d*d*d*b+28323/640*LAMBDA*LAMBDA*beta*d*beta*d*beta*d*d*b+107446159/1152*LAMBDA*LAMBDA*beta*d*beta*d*beta*d*beta*d*b+14013/320*LAMBDA*LAMBDA*beta*d*beta*d*beta*d*beta*d*beta*b+0*LAMBDA*LAMBDA*c*beta*d*beta*b*b*b*b+67360369/92160*LAMBDA*LAMBDA*beta*b*beta*b*beta*b*beta*b*beta*b*beta+47547/2560*LAMBDA*LAMBDA*beta*d*beta*b*beta*b*beta*b*beta*b*beta+28053/640*LAMBDA*LAMBDA*beta*d*beta*b*beta*d*beta*b*beta*b*beta+6561/160*LAMBDA*LAMBDA*beta*d*beta*b*beta*d*beta*d*beta*b*beta+243/20*LAMBDA*LAMBDA*beta*d*beta*d*beta*d*beta*d*beta*b*beta+0*LAMBDA*LAMBDA*d*d*d*d*d+0*LAMBDA*LAMBDA*beta*d*d*d*d*d+0*LAMBDA*LAMBDA*beta*d*beta*d*d*d*d+0*LAMBDA*LAMBDA*beta*d*beta*d*beta*d*d*d+0*LAMBDA*LAMBDA*beta*d*beta*d*beta*d*beta*d*d+0*LAMBDA*LAMBDA*beta*d*beta*d*beta*d*beta*d*beta*d+0*LAMBDA*LAMBDA*beta*d*beta*d*beta*d*beta*d*beta*d*beta+(-17637133/22544384)*LAMBDA*LAMBDA*b*b*b*b*b*b+12879/655360*LAMBDA*LAMBDA*b*beta*b*b*b*b*b+65853/655360*LAMBDA*LAMBDA*b*beta*b*b*beta*b*b*b+22113/81920*LAMBDA*LAMBDA*b*beta*b*b*beta*beta*b*b*b+4941/163840*LAMBDA*LAMBDA*b*d*b*b*b*b+6723/20480*LAMBDA*LAMBDA*b*beta*d*b*b*b*b+17323237/368640*LAMBDA*LAMBDA*b*beta*d*b*beta*b*b*b+529391/1440*LAMBDA*LAMBDA*b*beta*d*b*beta*beta*b*b*b+34147109/1474560*LAMBDA*LAMBDA*b*beta*b*b*beta*beta*b*b*beta*b+17213887/92160*LAMBDA*LAMBDA*b*beta*d*b*beta*beta*b*b*beta*b+12879/40960*LAMBDA*LAMBDA*b*beta*b*b*beta*beta*b*b*beta*beta*b+4293/1280*LAMBDA*LAMBDA*b*beta*d*b*beta*beta*b*b*beta*beta*b+8937/81920*LAMBDA*LAMBDA*b*d*b*d*b*b+23679/20480*LAMBDA*LAMBDA*b*beta*d*b*d*b*b+17711551/184320*LAMBDA*LAMBDA*b*beta*d*b*beta*d*b*b+17048161/23040*LAMBDA*LAMBDA*b*beta*d*b*beta*beta*d*b*b+17478271/46080*LAMBDA*LAMBDA*b*beta*d*b*beta*beta*d*b*beta*b+13311/1280*LAMBDA*LAMBDA*b*beta*d*b*beta*beta*d*b*beta*beta*b+4131/40960*LAMBDA*LAMBDA*b*beta*b*b*beta*beta*b*b*beta*beta*b*beta+2511/2560*LAMBDA*LAMBDA*b*beta*d*b*beta*beta*b*b*beta*beta*b*beta+3753/1280*LAMBDA*LAMBDA*b*beta*d*b*beta*beta*d*b*beta*beta*b*beta+297/2048*LAMBDA*LAMBDA*b*d*b*d*b*d+243/160*LAMBDA*LAMBDA*b*beta*d*b*d*b*d+16821/2560*LAMBDA*LAMBDA*b*beta*d*b*beta*d*b*d+2403/160*LAMBDA*LAMBDA*b*beta*d*b*beta*beta*d*b*d+33664511/5760*LAMBDA*LAMBDA*b*beta*d*b*beta*beta*d*b*beta*d+513/40*LAMBDA*LAMBDA*b*beta*d*b*beta*beta*d*b*beta*beta*d+567/160*LAMBDA*LAMBDA*b*beta*d*b*beta*beta*d*b*beta*beta*d*beta+135/2048*LAMBDA*LAMBDA*b*d*d*d*d*b+351/512*LAMBDA*LAMBDA*b*beta*d*d*d*d*b+7533/2560*LAMBDA*LAMBDA*b*beta*d*beta*d*d*d*b+2133/320*LAMBDA*LAMBDA*b*beta*d*beta*d*beta*d*d*b+5373/640*LAMBDA*LAMBDA*b*beta*d*beta*d*beta*d*beta*d*b+891/160*LAMBDA*LAMBDA*b*beta*d*beta*d*beta*d*beta*d*beta*b+0*LAMBDA*LAMBDA*d*d*d*d*d*b+0*LAMBDA*LAMBDA*d*beta*d*d*d*d*b+0*LAMBDA*LAMBDA*d*beta*d*beta*d*d*d*b+0*LAMBDA*LAMBDA*d*beta*d*beta*d*beta*d*d*b+0*LAMBDA*LAMBDA*d*beta*d*beta*d*beta*d*beta*d*b+0*LAMBDA*LAMBDA*d*beta*d*beta*d*beta*d*beta*d*beta*b+243/160*LAMBDA*LAMBDA*b*beta*d*beta*d*beta*d*beta*d*beta*b*beta+0*LAMBDA*LAMBDA*d*beta*d*beta*d*beta*d*beta*d*beta*b*beta+0*LAMBDA*LAMBDA*d*d*d*d*d*d+0*LAMBDA*LAMBDA*d*beta*d*d*d*d*d+0*LAMBDA*LAMBDA*d*beta*d*beta*d*d*d*d+0*LAMBDA*LAMBDA*d*beta*d*beta*d*beta*d*d*d+0*LAMBDA*LAMBDA*d*beta*d*beta*d*beta*d*beta*d*d+0*LAMBDA*LAMBDA*d*beta*d*beta*d*beta*d*beta*d*beta*d+0*LAMBDA*LAMBDA*d*beta*d*beta*d*beta*d*beta*d*beta*d*beta+0+(-524288/45)*LAMBDA*LAMBDA*beta*b*b*c+(-1048576/45)*LAMBDA*LAMBDA*beta*b*beta*beta*b*b*beta*c+(-131072/45)*LAMBDA*LAMBDA*beta*b*b*beta*beta*b*b*c+131072/45*LAMBDA*LAMBDA*beta*b*b*b*c");
//        t = CollectPowers.INSTANCE.transform(t);
        System.out.println(t);
    }
}
