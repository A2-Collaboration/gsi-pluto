{
    
    makeDistributionManager()->Startup("my = 0.2;eventnumber=0");

    PBatch * batch = new PBatch();

    batch->AddCommand("eventnumber=eventnumber+1");    
    batch->AddCommand("echo vorsub");
    batch->AddCommand("gosub sub");
    batch->AddCommand("echo nachsub");
    batch->AddCommand("xx: my = my + 0.1 ; echo the variable my is now:$my, eventnumber=$eventnumber ;  if (my < 0.5);  gosub sub; goto(xx)  ");

    batch->AddCommand("exit");
    //batch->AddCommand("goto end");
    batch->AddCommand("sub:");
    batch->AddCommand("echo sub $my");
    batch->AddCommand("gosub sub2");
    batch->AddCommand("return");
    
    batch->AddCommand("sub2:");
    batch->AddCommand("echo sub2");
    batch->AddCommand("return");

    batch->AddCommand("end:");
    batch->Execute();
    

}
