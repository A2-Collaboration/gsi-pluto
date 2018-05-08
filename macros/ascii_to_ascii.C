{

//Lets say that the file "multiple.lst" looks as follows:
//1,2,3
//3,4,5
//6,7,8
//9,1,2
//


    PReaction r;
    
    r.Input("multiple.lst", "readline{@a1,@a2,@a3};");
    r.Do("readline{@b1,@b2,@b3};"); 
    r.CloseFile();
    r.Do("c1=a1+b1; c2=a2+b2; c3=a3+b3;");
    r.Output("merged.lst","echo $c1,$c2,$c3");

    r->Print();

    cout << r.Loop() << " events converted" << endl;

}
