package cc.redberry.physics.OneLoopAction;

import cc.redberry.core.context.CC;
import cc.redberry.core.tensor.Expression;
import cc.redberry.core.tensor.SimpleTensor;
import cc.redberry.core.transformations.RenameConflictingIndices;
import cc.redberry.physics.util.SqrSubs;
import cc.redberry.transformation.CalculateNumbers;
import cc.redberry.transformation.ExpandBrackets;
import cc.redberry.transformation.Transformer;
import cc.redberry.transformation.collect.CollectFactory;
import cc.redberry.transformation.collect.CollectPowers;
import cc.redberry.transformation.concurrent.EACScalars;
import cc.redberry.transformation.contractions.IndicesContractionsTransformation;

public class OneLoop extends MainTensors {
    public static Expression Gamma;
    public static Expression Lambda;

    public static void main(String[] args) {
        K_2 = new Expression("K^{\\mu\\nu}_\\alpha^\\beta=g^{\\mu\\nu}*d_{\\alpha}^{\\beta}-\\lambda/2*(g^{\\mu\\beta}*d_\\alpha^\\nu+g^{\\nu\\beta}*d_\\alpha^\\mu)");
        KINV = new Expression("KINV_\\alpha^\\beta=d_\\alpha^\\beta-\\gamma*n_\\alpha*n^\\beta");
        W = new Expression("W^{\\alpha}_{\\beta}=P^{\\alpha}_{\\beta}+\\lambda/2*R^\\alpha_\\beta");
        calculate(K_2, KINV, W);
    }

    public OneLoop() {
         start();
    }

    public static OneLoop calculate(Expression K_2_in, Expression KINV_in, Expression W_in) {
        L = new Expression("L=2");
        Gamma = new Expression("\\gamma=\\lambda/(1-\\lambda)");
        Lambda = new Expression("\\lambda=0");
        K_2 = K_2_in;
        KINV = KINV_in;
        W = W_in;
        HATS_0 = new Expression("HATS=0");
        HATS_1 = new Expression("HATS_\\alpha=0");
        HATS_2 = new Expression("HATS_{\\alpha\\beta}=0");
        HATS_3 = new Expression("HATS_{\\alpha\\beta\\gamma}=0");
        HATS_4 = new Expression("HATS_{\\alpha\\beta\\gamma\\nu}=0");
        NABLAS_0 = new Expression("NABLAS=0");
        NABLAS_1 = new Expression("NABLAS_\\alpha=0");
        NABLAS_2 = new Expression("NABLAS_{\\alpha\\beta}=0");
        NABLAS_3 = new Expression("NABLAS_{\\alpha\\beta\\gamma}=0");
        NABLAS_4 = new Expression("NABLAS_{\\alpha\\beta\\gamma\\nu}=0");
        HATN_0 = new Expression("HATN=0");
        HATN_1 = new Expression("HATN_\\alpha=0");
        HATN_2 = new Expression("HATN_{\\alpha\\beta}=0");
        HATN_3 = new Expression("HATN_{\\alpha\\beta\\gamma}=0");
        HATN_4 = new Expression("HATN_{\\alpha\\beta\\gamma\\nu}=0");
        HATM_0 = new Expression("HATM=0");
        HATM_1 = new Expression("HATM_\\alpha=0");
        HATM_2 = new Expression("HATM_{\\alpha\\beta}=0");
        HATM_3 = new Expression("HATM_{\\alpha\\beta\\gamma}=0");
        HATM_4 = new Expression("HATM_{\\alpha\\beta\\gamma\\nu}=0");

        INPUTs = new Expression[]{K_2, KINV, W, HATS_0, HATS_1, HATS_2, HATS_3, HATS_4, HATN_0, HATN_1, HATN_2, HATN_3, HATN_4, HATM_0, HATM_1, HATM_2, HATM_3, HATM_4, NABLAS_0, NABLAS_1, NABLAS_2, NABLAS_3, NABLAS_4};

        return new OneLoop();
    }

    public void start() {
        System.out.println("Calculating...");
        Gamma.eval(Lambda);
        K_2.eval(Lambda);
        KINV.eval(Gamma);
        W.eval(Lambda);
        for (Expression ex : ALL)
            ex.eval(L, CalculateNumbers.INSTANCE);
        for (Expression ex : INPUTs)
            ex.eval(L, CalculateNumbers.INSTANCE);
        evalHats();
        //System.out.println("---------DALTAs-----------");
        evalDeltas();
        //System.out.println("---------TERMs-----------");
        evalTerms();
       // System.out.println("----------ACTION----------");
        evalAction();
        //latex.push(ACTION.toString());
        //latex.output();
    }

    public final void evalHats() {
        for (Expression hat : HATs) {
            for (Expression in : INPUTs)
                hat.eval(in);
            hat.eval(
                    new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    KRONECKER_DIMENSION,
                    CollectFactory.createCollectEqualTerms(),
                    CalculateNumbers.INSTANCE,
                    EACScalars.getTransformer(),
                    CalculateNumbers.INSTANCE);
        }
    }

    private void evalDeltas() {
        for (Expression delta : DELTAs) {
            for (Expression hat : HATs)
                delta.eval(hat);
            delta.eval(
                    new Transformer(new POO()),
                    CalculateNumbers.INSTANCE,
                    new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    KRONECKER_DIMENSION,
                    CollectFactory.createCollectAllScalars(),
                    CalculateNumbers.INSTANCE);
        }
    }

    public final void evalTerms() {
        for (Expression e : TERMs) {
            for (Expression d : DELTAs)
                e.eval(
                        d.asSubstitution());
            for (Expression d : HATs)
                e.eval(
                        d.asSubstitution());

            for (Expression d : INPUTs)
                e.eval(
                        d.asSubstitution(),
                        CalculateNumbers.INSTANCE);
            e.eval(
                    new Transformer(new POO()),
                    CalculateNumbers.INSTANCE,
                    new Transformer(RenameConflictingIndices.INSTANCE),
                    new Transformer(ExpandBrackets.EXPAND_ALL),
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    new SqrSubs((SimpleTensor) CC.parse("n_{\\alpha}")),
                    KRONECKER_DIMENSION,
                    CalculateNumbers.INSTANCE,
                    new Transformer(CollectPowers.INSTANCE),
                    CollectFactory.createCollectEqualTerms(),
                    CalculateNumbers.INSTANCE,
                    Averaging.INSTANCE,
                    new Transformer(ExpandBrackets.EXPAND_ALL),
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    KRONECKER_DIMENSION,
                    CalculateNumbers.INSTANCE,
                    new Expression("R_{\\mu \\nu}^{\\mu}_{\\alpha} = R_{\\nu\\alpha}"),
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    new Expression("R_{\\mu}^{\\mu} = R"),
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    CollectFactory.createCollectAllEqualTerms(),
                    new Transformer(new POO()),
                    CalculateNumbers.INSTANCE,
                    new Expression("R_{\\mu\\nu}^{\\alpha}_{\\alpha}=0"),
                    CalculateNumbers.INSTANCE,
                    CollectFactory.createCollectAllEqualTerms(),
                    CalculateNumbers.INSTANCE,
                    new Expression("F_\\mu^\\mu_\\alpha\\beta=0"),
                    new Expression("F_\\alpha\\beta_\\mu^\\mu=0"),
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    new Expression("R_{\\mu\\nu\\alpha\\beta}*R^{\\mu\\alpha\\nu\\beta}=(1/2)*R_{\\mu\\nu\\alpha\\beta}*R^{\\mu\\nu\\alpha\\beta}"),
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    new Expression("R_{\\mu\\nu\\alpha\\beta}*R^{\\mu\\nu\\alpha\\beta}=4*R_{\\mu\\nu}*R^{\\mu\\nu}-R*R"),
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    new Expression("R_{\\mu \\nu}^{\\mu}_{\\alpha} = R_{\\nu\\alpha}"),
                    new Transformer(ExpandBrackets.EXPAND_ALL),
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    new Expression("R_{\\mu\\nu}^{\\alpha}_{\\alpha}=0"),
                    CollectFactory.createCollectAllEqualTerms(),
                    CalculateNumbers.INSTANCE
            );
        }

    }

    public final void evalAction() {

        ACTION.eval(
                Flat, WR, SR, SSR, FF, FR, RR);
        ACTION.eval(
                new Expression("F_\\mu\\nu\\alpha\\beta=R_\\mu\\nu\\alpha\\beta"),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                new Expression("R_{\\mu\\nu\\alpha\\beta}*R^{\\mu\\alpha\\nu\\beta}=(1/2)*R_{\\mu\\nu\\alpha\\beta}*R^{\\mu\\nu\\alpha\\beta}"),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                new Expression("R_{\\mu\\nu\\alpha\\beta}*R^{\\mu\\nu\\alpha\\beta}=4*R_{\\mu\\nu}*R^{\\mu\\nu}-R*R"),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                new Expression("R_{\\mu \\nu}^{\\mu}_{\\alpha} = R_{\\nu\\alpha}"),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                new Expression("R_{\\mu\\nu}^{\\alpha}_{\\alpha}=0"),
                new Transformer(ExpandBrackets.EXPAND_ALL),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                CollectFactory.createCollectAllEqualTerms(),
                new Transformer(CollectFactory.createCollectAllScalars()),
                CalculateNumbers.INSTANCE);
    }

}
