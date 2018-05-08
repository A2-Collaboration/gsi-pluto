




{

    makeDistributionManager()->Exec("elementary");
    makeDistributionManager()->LinkDB();
    makeDistributionManager()->GetDistribution("D+_prod_tcross")->SetRange(2.,5.);
    makeDistributionManager()->GetDistribution("D+_prod_tcross")->Draw();


}
