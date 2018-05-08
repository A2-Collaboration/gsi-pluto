{

    Int_t nEvents;
    nEvents = 1000;

    

    //this must come in the very beginning, because here we "freeze-out"
    //the physics
    PDecayManager *p_p = new PDecayManager;

    

    listParticle("Lambda1405");
    listParticle("Sigma1385+");
    listParticle("Sigma1385-");
    listParticle("Sigma13850");

    


    

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //calculation of the branching ratio

    double channel_1,channel_2, channel_3,channel_4,channel_5,channel_6,channel_7,channel_8,channel_9,channel_10,channel_11,channel_12,channel_13;
    double channel_14,channel_16,channel_17,channel_18,channel_19,channel_22,channel_23,channel_24,channel_25,channel_26,channel_27;
    double channel_28,channel_29,channel_30,channel_31,channel_32,channel_33,channel_34,channel_35,channel_36,channel_37,channel_38;
    double channel_39,channel_41,channel_42,channel_43,channel_46,channel_47,channel_48,channel_50,channel_51,channel_53,channel_54;
    double channel_55,channel_56,channel_57,channel_59,channel_62,channel_63,channel_64,channel_65,channel_66,channel_68,channel_73;
    double channel_75,channel_76,channel_77,channel_78,channel_81,channel_83,channel_84,channel_85,channel_86,channel_87,channel_91;
    double channel_92,channel_93,channel_96,channel_97,channel_99,channel_100,channel_101,channel_102,channel_103,channel_104;
    double channel_105,channel_106,channel_107,channel_108,channel_109;



    double total;
    double total_sum;

    double width_channel_1,width_channel_2,width_channel_3,width_channel_4,width_channel_5,width_channel_6,width_channel_7,width_channel_8,width_channel_9;
    double width_channel_10,width_channel_11,width_channel_12,width_channel_13,width_channel_14,width_channel_16,width_channel_17;
    double width_channel_18,width_channel_19,width_channel_22,width_channel_23,width_channel_24,width_channel_25,width_channel_26;
    double width_channel_27,width_channel_28,width_channel_29,width_channel_30,width_channel_31,width_channel_32,width_channel_33;
    double width_channel_34,width_channel_35,width_channel_36,width_channel_37,width_channel_38,width_channel_39,width_channel_41;
    double width_channel_42,width_channel_43,width_channel_46,width_channel_47,width_channel_48,width_channel_50,width_channel_51;
    double width_channel_53,width_channel_54,width_channel_55,width_channel_56,width_channel_57,width_channel_59,width_channel_62;
    double width_channel_63,width_channel_64,width_channel_65,width_channel_66,width_channel_68,width_channel_73,width_channel_75;
    double width_channel_76,width_channel_77,width_channel_78,width_channel_81,width_channel_83,width_channel_84,width_channel_85;
    double width_channel_86,width_channel_87,width_channel_91,width_channel_92,width_channel_93,width_channel_96,width_channel_97;
    double width_channel_99,width_channel_100,width_channel_101,width_channel_102,width_channel_103,width_channel_104,width_channel_105;
    double width_channel_106,width_channel_107,width_channel_108,width_channel_109;


    channel_1=0.012446;    //neuer Deuteron Kanal
    channel_2=8.75415;
    channel_3=4.939;         //4.89153;
    channel_4=3.125;         //4.11721;
    channel_5=3.23167;       //3.50459;
    channel_6=2.80371;
    channel_7=0.0252;        //2.50808; eta'
    channel_8=1.86992;
    channel_9=1.82383;
    channel_10=1.28704;
    channel_11=1.606;         //1.15488;
    channel_12=0.15271;       //#1.2543#;      //0.96555; **
    channel_13=1.14247;      //0.778186;
    channel_14=0.83239;       //0.745355;
    channel_16=0.607986;
    channel_17=0.559463;
    channel_18=0.446415;
    channel_19=0.342269;
    channel_22=0.153236;
    channel_23=0.320938;
    channel_24=0.256578;
    channel_25=0.00010838;     //#0.0152#     //0.162354; **
    channel_26=0.59145;       //bei selektiert -->0.61554     //0.144175; ******* sieht scheiﬂe aus
    channel_27=0.135824;
    channel_28=0.130073;
    channel_29=0.116638;
    channel_30=0.705;        //0.115416;
    channel_31=0.27209;       //0.112524;
    channel_32=0.111833;
    channel_33=0.0911551;
    channel_34=0.0811768;
    channel_35=0.0810533;


    channel_36=0.0785387;
    channel_37=0.05609;
    channel_38=0.0523881;
    channel_39=0.0513331;
    channel_41=0.0461778;
    channel_42=0.0443434;
    channel_43=0.04673;        //0.0433381;   //0.05  benjamin sailer
    channel_46=0.0412845;
    channel_47=0.0390705;
    channel_48=0.0372832;
    channel_50=0.0369692;
    channel_51=0.0341877;
    channel_53=0.0329389;
    channel_54=0.0286673;
    channel_55=0.0276463;
    channel_56=0.0256681;
    channel_57=0.0239695;
    channel_59=0.0215231;
    channel_62=0.0172551;
    channel_63=0.0171834;
    channel_64=0.0155394;
    channel_65=0.014081;
    channel_66=0.0123809;
    channel_68=0.0103284;
    channel_73=0.0073353;
    channel_75=0.00604724;
    //channel_76=0.00585238;
    channel_77=0.00534218;
    channel_78=0.00532157;
    channel_81=0.00508021;
    //channel_83=0.00459111;


    channel_84=0.00453141;
    channel_85=0.00446203;
    channel_86=0.00406263;
    channel_87=0.00375049;
    channel_91=0.00227592;
    channel_92=0.00225889;
    channel_93=0.00202397;
    channel_96=0.00191676;
    //channel_97=0.00179125;
    channel_99=0.00132698;
    channel_100=0.00132376;
    channel_101=0.00116988;
    channel_102=0.00112839;
    channel_103=0.00102256;
    //channel_104=0.000843518;
    channel_105=0.000583393;
    channel_106=0.000372088;
    channel_107=0.000354397;
    //channel_108=0.000165741;
    channel_109=0.000039595;



    //total = channel_1 + channel_2;

    total =   channel_2+channel_3+channel_4+channel_5+channel_6+channel_7+channel_8+channel_9+channel_10+channel_11+channel_12;
    total +=  channel_13+channel_14+channel_16+channel_17+channel_18+channel_19+channel_22+channel_23+channel_24+channel_25;
    total +=  channel_26+channel_27+channel_28+channel_29+channel_30+channel_31+channel_32+channel_33+channel_34+channel_35;

    total +=  channel_1+channel_36+channel_37+channel_38+channel_39+channel_41+channel_42+channel_43+channel_46+channel_47;
    total +=  channel_48+channel_50+channel_51+channel_53+channel_54+channel_55+channel_56+channel_57+channel_59+channel_62+channel_63;
    total +=  channel_64+channel_65+channel_66+channel_68+channel_73+channel_75+channel_77+channel_78+channel_81;


    total +=  channel_84+channel_85+channel_86+channel_87+channel_91+channel_92+channel_93+channel_96;
    total +=  channel_99+channel_100+channel_101+channel_102+channel_103+channel_105+channel_106+channel_107;
    total +=  channel_109;


    //     ++channel_76+channel_108+channel_83+channel_97+channel_104
    //PReaction: insufficient energy
    //retval: 6


    cout<<"total :"<<total<<endl;

    width_channel_1 = channel_1/total;
    width_channel_2 = channel_2/total;
    width_channel_3 = channel_3/total;
    width_channel_4 = channel_4/total;
    width_channel_5 = channel_5/total;
    width_channel_6 = channel_6/total;
    width_channel_7 = channel_7/total;
    width_channel_8 = channel_8/total;
    width_channel_9 = channel_9/total;
    width_channel_10 = channel_10/total;
    width_channel_11 = channel_11/total;
    width_channel_12 = channel_12/total;
    width_channel_13 = channel_13/total;
    width_channel_14 = channel_14/total;
    width_channel_16 = channel_16/total;
    width_channel_17 = channel_17/total;
    width_channel_18 = channel_18/total;
    width_channel_19 = channel_19/total;
    width_channel_22 = channel_22/total;
    width_channel_23 = channel_23/total;
    width_channel_24 = channel_24/total;
    width_channel_25 = channel_25/total;
    width_channel_26 = channel_26/total;
    width_channel_27 = channel_27/total;
    width_channel_28 = channel_28/total;
    width_channel_29 = channel_29/total;
    width_channel_30 = channel_30/total;
    width_channel_31 = channel_31/total;
    width_channel_32 = channel_32/total;
    width_channel_33 = channel_33/total;
    width_channel_34 = channel_34/total;
    width_channel_35 = channel_35/total;
    width_channel_36 = channel_36/total;
    width_channel_37 = channel_37/total;
    width_channel_38 = channel_38/total;
    width_channel_39 = channel_39/total;
    width_channel_41 = channel_41/total;
    width_channel_42 = channel_42/total;
    width_channel_43 = channel_43/total;
    width_channel_46 = channel_46/total;
    width_channel_47 = channel_47/total;
    width_channel_48 = channel_48/total;
    width_channel_50 = channel_50/total;
    width_channel_51 = channel_51/total;
    width_channel_53 = channel_53/total;
    width_channel_54 = channel_54/total;
    width_channel_55 = channel_55/total;
    width_channel_56 = channel_56/total;
    width_channel_57 = channel_57/total;
    width_channel_59 = channel_59/total;
    width_channel_62 = channel_62/total;
    width_channel_63 = channel_63/total;
    width_channel_64 = channel_64/total;
    width_channel_65 = channel_65/total;
    width_channel_66 = channel_66/total;
    width_channel_68 = channel_68/total;
    width_channel_73 = channel_73/total;
    width_channel_75 = channel_75/total;
    //width_channel_76 = channel_76/total;
    width_channel_77 = channel_77/total;
    width_channel_78 = channel_78/total;
    width_channel_81 = channel_81/total;
    //width_channel_83 = channel_83/total;
    width_channel_84 = channel_84/total;
    width_channel_85 = channel_85/total;
    width_channel_86 = channel_86/total;
    width_channel_87 = channel_87/total;
    width_channel_91 = channel_91/total;
    width_channel_92 = channel_92/total;
    width_channel_93 = channel_93/total;
    width_channel_96 = channel_96/total;
    //width_channel_97 = channel_97/total;
    width_channel_99 = channel_99/total;
    width_channel_100 = channel_100/total;
    width_channel_101 = channel_101/total;
    width_channel_102 = channel_102/total;
    width_channel_103 = channel_103/total;
    //width_channel_104 = channel_104/total;
    width_channel_105 = channel_105/total;
    width_channel_106 = channel_106/total;
    width_channel_107 = channel_107/total;
    //width_channel_108 = channel_108/total;
    width_channel_109 = channel_109/total;



    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Float_t Eb    = 3.5;         // beam energy in AGeV

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Some Particles for the reaction

    PParticle *p_b=new PParticle("p",Eb);  //Beam proton
    PParticle *p_t=new PParticle("p");     //Target proton
    PParticle *q=new PParticle(*p_b+*p_t);



    PParticle *Lambda1405=new PParticle("Lambda1405");//Lambda1405
    PParticle *Sigma1385P=new PParticle("Sigma1385+");//Lambda1405
    PParticle *Sigma1385M=new PParticle("Sigma1385-");//Lambda1405
    PParticle *Sigma1385N=new PParticle("Sigma13850");//Lambda1405
    PParticle *Lambda=new PParticle("Lambda");        //Lambda
    PParticle *SigmaP=new PParticle("Sigma+");        //Lambda
    PParticle *SigmaM=new PParticle("Sigma-");        //Lambda
    PParticle *SigmaN=new PParticle("Sigma0");        //Lambda
    PParticle *kplus=new PParticle("K+");             //K+ meson
    PParticle *kminus=new PParticle("K-");            //K- meson
    PParticle *K0S=new PParticle("K0S");              //K0 short
    PParticle *K0S2=new PParticle("K0S");              //K0 short
    PParticle *Proton1=new PParticle("p");            //Proton 1
    PParticle *Proton2=new PParticle("p");            //Proton 2
    PParticle *Neutron=new PParticle("n");            //Neutron
    PParticle *PiPlus1=new PParticle("pi+");          //K0 short
    PParticle *PiPlus2=new PParticle("pi+");
    PParticle *PiPlus3=new PParticle("pi+");
    PParticle *PiPlus4=new PParticle("pi+");
    PParticle *PiPlus5=new PParticle("pi+");
    PParticle *PiMinus1=new PParticle("pi-");           //PiMinus
    PParticle *PiMinus2=new PParticle("pi-");           //PiMinus
    PParticle *PiMinus3=new PParticle("pi-");           //PiMinus
    PParticle *PiMinus4=new PParticle("pi-");           //PiMinus
    PParticle *PiNull1=new PParticle("pi0");            //Neutro
    PParticle *PiNull2=new PParticle("pi0");            //Neutro
    PParticle *DP=new PParticle("D+");
    PParticle *DPP=new PParticle("D++");
    PParticle *DM=new PParticle("D-");
    PParticle *DNull=new PParticle("D0");
    PParticle *etaP=new PParticle("eta'");
    PParticle *eta=new PParticle("eta");
    PParticle *omega=new PParticle("w");
    PParticle *rhoM=new PParticle("rho-");
    PParticle *rhoP=new PParticle("rho+");
    PParticle *rhoNull=new PParticle("rho0");
    PParticle *N1440Null=new PParticle("NP110");
    PParticle *N1520Null=new PParticle("ND130");
    PParticle *N1535Null=new PParticle("NS110");
    PParticle *N1440P=new PParticle("NP11+");
    PParticle *N1520P=new PParticle("ND13+");
    PParticle *N1535P=new PParticle("NS11+");
    PParticle *d=new PParticle("d");
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Set decay manager



    PDecayChannel *c = new PDecayChannel;



     p_p->SetDefault("Lambda1405");
     p_p->SetDefault("Sigma1385+");
     p_p->SetDefault("Sigma1385-");
     p_p->SetDefault("Sigma13850");
     p_p->SetDefault("D++");
     p_p->SetDefault("D+");
     p_p->SetDefault("D-");
     p_p->SetDefault("D0");
     p_p->SetDefault("rho-");
     p_p->SetDefault("rho+");
     p_p->SetDefault("rho0");
     p_p->SetDefault("NP110");
     p_p->SetDefault("ND130");
     p_p->SetDefault("NS110");
     p_p->SetDefault("NP11+");
     p_p->SetDefault("ND13+");
     p_p->SetDefault("NS11+");
     p_p->SetDefault("w");
     p_p->SetDefault("eta'");
     p_p->SetDefault("dimuon");
     p_p->SetDefault("dilepton");

     p_p->SetDefault("eta");
     //p_p->SetDefault("pi0");

    //Add all participated channels

    PParticle * channel0[] ={d,PiPlus1};
    c->AddChannel(width_channel_1,2,channel0);
    //c->AddChannel(width_channel_1,"d","pi+");

#if 1
    PParticle * channel1[] ={Proton1,Neutron,PiPlus1};
    c->AddChannel(width_channel_2,3,channel1);                       // include decay modes

    
    PParticle * channel2[] ={Neutron,DPP};
    c->AddChannel(width_channel_3,2,channel2);

    PParticle * channel3[] ={PiPlus1,PiNull1,Proton1,Neutron};
    c->AddChannel(width_channel_4,4,channel3);                       // include decay modes


    PParticle * channel4[] ={PiNull1,Proton1,Proton2};
    c->AddChannel(width_channel_5,3,channel4);                       // include decay modes
    PParticle * channel5[] ={PiPlus1,PiMinus1,Proton1,Proton2};
    c->AddChannel(width_channel_6,4,channel5);                       // include decay modes
    PParticle * channel6[] ={etaP,Proton1,Proton2};
    c->AddChannel(width_channel_7,3,channel6);                      // include decay modes
    PParticle * channel7[] ={Proton1,Neutron,PiNull1,PiNull2,PiPlus1};
    c->AddChannel(width_channel_8,5,channel7);                      // include decay modes
    PParticle * channel8[] ={Proton1,Neutron,PiMinus1,PiPlus2,PiPlus1};
    c->AddChannel(width_channel_9,5,channel8);                      // include decay modes
    PParticle * channel9[] ={Proton1,Proton2,PiMinus1,PiPlus1,PiNull1};
    c->AddChannel(width_channel_10,5,channel9);                      // include decay modes

    PParticle * channel10[] ={Proton1,DP};
    c->AddChannel(width_channel_11,2,channel10);                     // include decay modes
    PParticle * channel11[] ={Proton1,Proton2,PiNull1,PiNull2};
    c->AddChannel(width_channel_12,4,channel11);                     // include decay modes
    PParticle * channel12[] ={Neutron,DPP,PiMinus1,PiPlus2};
    c->AddChannel(width_channel_13,4,channel12);                     // include decay modes
    PParticle * channel13[] ={DM,Proton1,PiPlus1,PiPlus2};
    c->AddChannel(width_channel_14,4,channel13);                     // include decay modes

    PParticle * channel14[] ={Proton1,N1520P};
    c->AddChannel(width_channel_16,2,channel14);                     // include decay modes

    PParticle * channel15[] ={DPP,DM,PiPlus1};
    c->AddChannel(width_channel_17,3,channel15);

    PParticle * channel16[] ={Proton1,N1535P};
    c->AddChannel(width_channel_18,2,channel16);                    // include decay modes


    PParticle * channel17[] ={Proton1,PiNull1,N1520P};                //???????????????????
    c->AddChannel(width_channel_19,3,channel17);                     // include decay modes


    PParticle * channel18[] ={Proton1,PiNull1,N1440P};
    c->AddChannel(width_channel_22,3,channel18);                    // include decay modes

    PParticle * channel19[] ={Proton1,PiMinus1,DPP};
    c->AddChannel(width_channel_23,3,channel19);                    // include decay modes

    PParticle * channel20[] ={Proton1,PiMinus1,DPP,PiMinus2,PiPlus1};
    c->AddChannel(width_channel_24,5,channel20);                    // include decay modes
    PParticle * channel21[] ={d,PiPlus1,PiNull1};
    c->AddChannel(width_channel_25,3,channel21);                     // include decay modes
    PParticle * channel22[] ={Proton1,N1440P};
    c->AddChannel(width_channel_26,2,channel22);                     // include decay modes
    PParticle * channel23[] ={Proton1,Proton2,PiPlus1,PiMinus1,PiPlus2,PiMinus2};
    c->AddChannel(width_channel_27,6,channel23);                     // include decay modes
    PParticle * channel24[] ={Proton1,PiNull1,N1535P};
    c->AddChannel(width_channel_28,3,channel24);                     // include decay modes
    PParticle * channel25[] ={Proton1,Proton2,PiPlus1,rhoM};
    c->AddChannel(width_channel_29,4,channel25);                       // include decay modes
    PParticle * channel26[] ={Proton1,DPP,PiNull1,PiMinus1};
    c->AddChannel(width_channel_30,4,channel26);                     // include decay modes
    PParticle * channel27[] ={Proton1,Proton2,omega};
    c->AddChannel(width_channel_31,3,channel27);                     // include decay modes
    PParticle * channel28[] ={Neutron,DPP,PiPlus1,PiMinus1,PiPlus2,PiMinus2};
    c->AddChannel(width_channel_32,6,channel28);                     // include decay modes

    PParticle * channel29[] ={Proton1,DPP,PiPlus1,PiNull1,PiMinus1,PiMinus2};
    c->AddChannel(width_channel_33,6,channel29);                     // include decay modes
    PParticle * channel30[] ={Proton1,N1520Null,PiPlus1,PiNull1};
    c->AddChannel(width_channel_34,4,channel30);                     // include decay modes
    PParticle * channel31[] ={Proton1,Proton2,PiMinus1,rhoP};
    c->AddChannel(width_channel_35,4,channel31);                     // include decay modes




    PParticle * channel32[] ={Proton1,Proton2,PiNull1,rhoNull};
    c->AddChannel(width_channel_36,4,channel32);                     // include decay modes
    PParticle * channel33[] ={Proton1,DM,PiPlus1,PiPlus2,PiPlus3,PiMinus1};
    c->AddChannel(width_channel_37,6,channel33);                     // include decay modes
    PParticle * channel34[] ={Neutron,DPP,rhoNull};
    c->AddChannel(width_channel_38,3,channel34);                     // include decay modes
    PParticle * channel35[] ={PiPlus1,d,PiPlus2,PiNull1,PiMinus1};
    c->AddChannel(width_channel_39,5,channel35);
    PParticle * channel36[] ={Proton1,Proton2,PiPlus1,PiPlus2,PiNull1,PiMinus1,PiMinus2};
    c->AddChannel(width_channel_41,7,channel36);                     // include decay modes
    PParticle * channel37[] ={Neutron,DP,rhoP};
    c->AddChannel(width_channel_42,3,channel37);                     // include decay modes
    PParticle * channel38[] ={Proton1,Proton2,eta};
    c->AddChannel(width_channel_43,3,channel38);                     // include decay modes
    PParticle * channel39[] ={d,PiPlus1,PiPlus2,PiMinus1};
    c->AddChannel(width_channel_46,4,channel39);                     // include decay modes

    PParticle * channel40[] ={Proton1,Proton2,PiPlus1,PiMinus1,rhoNull};
    c->AddChannel(width_channel_47,5,channel40);                     // include decay modes

    PParticle * channel41[] ={Proton1,Neutron,PiPlus1,PiPlus2,PiPlus3,PiMinus1,PiMinus2};
    c->AddChannel(width_channel_48,7,channel41);                     // include decay modes
    PParticle * channel42[] ={Proton1,Proton2,rhoNull};
    c->AddChannel(width_channel_50,3,channel42);                     // include decay modes
    PParticle * channel43[] ={Proton1,DPP,PiPlus1,PiPlus2,PiMinus1,PiMinus2,PiMinus3};
    c->AddChannel(width_channel_51,7,channel43);                     // include decay modes
    PParticle * channel44[] ={Proton1,N1520P,PiPlus1,PiMinus1};
    c->AddChannel(width_channel_53,4,channel44);                     // include decay modes
    PParticle * channel45[] ={Proton1,DPP,rhoM};
    c->AddChannel(width_channel_54,3,channel45);                     // include decay modes
    PParticle * channel46[] ={Lambda,Proton1,PiPlus1,K0S};
    c->AddChannel(width_channel_55,4,channel46);                     // include decay modes
    PParticle * channel47[] ={Proton1,DNull,rhoP};
    c->AddChannel(width_channel_56,3,channel47);                     // include decay modes
    PParticle * channel48[] ={SigmaP,Proton1,PiNull1};
    c->AddChannel(width_channel_57,3,channel48);                     // include decay modes
    PParticle * channel49[] ={Proton1,N1520P,PiPlus1,PiNull1,PiMinus1};
    c->AddChannel(width_channel_59,5,channel49);                     // include decay modes
    PParticle * channel50[] ={Proton1,Proton2,PiPlus1,PiNull1,PiMinus1,rhoNull};
    c->AddChannel(width_channel_62,6,channel50);                     // include decay modes
    PParticle * channel51[] ={PiPlus1,DPP,DM,PiPlus2,PiMinus1};
    c->AddChannel(width_channel_63,5,channel51);                     // include decay modes

    PParticle * channel52[] ={Proton1,Proton2,PiPlus1,omega,PiMinus1};
    c->AddChannel(width_channel_64,5,channel52);                     // include decay modes

    PParticle * channel53[] ={Proton1,N1520Null,PiPlus1,PiPlus2,PiMinus1};
    c->AddChannel(width_channel_65,5,channel53);                     // include decay modes
    PParticle * channel54[] ={SigmaN,Proton1,PiPlus1,K0S};
    c->AddChannel(width_channel_66,4,channel54);                     // include decay modes
    PParticle * channel55[] ={Lambda,DPP,K0S};
    c->AddChannel(width_channel_68,3,channel55);                     // include decay modes

    PParticle * channel56[] ={Proton1,Proton2,PiNull1,K0S,K0S2};
    c->AddChannel(width_channel_73,5,channel56);                     // include decay modes
    PParticle * channel57[] ={Proton1,Proton2,PiPlus1,PiPlus2,PiMinus1,rhoM};
    c->AddChannel(width_channel_75,6,channel57);                     // include decay modes
    //PParticle * channel58[] ={Proton1,Proton2,PiPlus1,PiPlus2,PiPlus3,PiNull1,PiMinus1,PiMinus2,PiMinus3};
    //c->AddChannel(width_channel_76,9,channel58);                     // include decay modes
    PParticle * channel59[] ={Sigma1385P,Proton1,K0S};
    c->AddChannel(width_channel_77,3,channel59);                     // include decay modes
    PParticle * channel60[] ={Proton1,Proton2,PiPlus1,PiMinus1,PiMinus2,rhoP};
    c->AddChannel(width_channel_78,6,channel60);                     // include decay modes
    PParticle * channel61[] ={Lambda,Neutron,PiPlus1,PiPlus2,K0S};
    c->AddChannel(width_channel_81,5,channel61);                     // include decay modes






    //PParticle * channel62[] ={Proton1,Neutron,PiPlus1,PiPlus2,PiPlus3,PiPlus4,PiMinus1,PiMinus2,PiMinus3};
    //c->AddChannel(width_channel_83,9,channel62);                     // include decay modes
    PParticle * channel63[] ={SigmaP,Neutron,PiPlus1,K0S};
    c->AddChannel(width_channel_84,4,channel63);                     // include decay modes

    PParticle * channel64[] ={Lambda,Proton1,PiPlus1,PiNull1,K0S};
    c->AddChannel(width_channel_85,5,channel64);                     // include decay modes
    PParticle * channel65[] ={SigmaP,Proton1,PiNull1,K0S};
    c->AddChannel(width_channel_86,4,channel65);                     // include decay modes
    PParticle * channel66[] ={Proton1,SigmaM,PiPlus1,PiPlus2,K0S};
    c->AddChannel(width_channel_87,5,channel66);                     // include decay modes
    PParticle * channel67[] ={Proton1,DPP,PiMinus1,omega};
    c->AddChannel(width_channel_91,4,channel67);                     // include decay modes

    PParticle * channel68[] ={Proton1,SigmaP,PiPlus1,PiMinus1,K0S};
    c->AddChannel(width_channel_92,5,channel68);                     // include decay modes
    PParticle * channel69[] ={Proton1,Proton2,PiPlus1,K0S,kminus};
    c->AddChannel(width_channel_93,5,channel69);                     // include decay modes
    PParticle * channel70[] ={Proton1,Neutron,PiPlus1,K0S,K0S2};
    c->AddChannel(width_channel_96,5,channel70);                     // include decay modes

    //PParticle * channel71[] ={Proton1,Proton2,PiPlus1,PiPlus2,PiPlus3,PiMinus1,PiMinus2,PiMinus3};
    //c->AddChannel(width_channel_97,8,channel71);                     // include decay modes
    PParticle * channel72[] ={Proton1,Proton2,K0S,K0S2};                                               // ???????????????????????
    c->AddChannel(width_channel_99,4,channel72);                     // include decay modes
    PParticle * channel73[] ={Lambda,Proton1,PiPlus1,PiPlus2,PiMinus1,K0S};
    c->AddChannel(width_channel_100,6,channel73);                     // include decay modes
    PParticle * channel74[] ={Proton1,SigmaM,PiPlus1,PiPlus2,PiNull1,K0S};
    c->AddChannel(width_channel_101,6,channel74);                     // include decay modes
    PParticle * channel75[] ={Lambda,Proton1,PiPlus1,PiPlus2,PiNull1,PiMinus1,K0S};
    c->AddChannel(width_channel_102,7,channel75);                     // include decay modes
    PParticle * channel76[] ={Proton1,Proton2,K0S,K0S2};
    c->AddChannel(width_channel_103,4,channel76);                     // include decay modes
    //PParticle * channel77[] ={Proton1,Neutron,PiPlus1,PiPlus2,PiPlus3,PiPlus4,PiPlus5,PiMinus1,PiMinus2,PiMinus3,PiMinus4};
    //c->AddChannel(width_channel_104,11,channel77);                    // include decay modes
    PParticle * channel78[] ={Neutron,SigmaM,PiPlus1,PiPlus2,PiPlus3,K0S};
    c->AddChannel(width_channel_105,6,channel78);                     // include decay modes
    PParticle * channel79[] ={Lambda,Neutron,PiPlus1,PiPlus2,PiPlus3,PiMinus1,K0S};
    c->AddChannel(width_channel_106,7,channel79);                     // include decay modes
    PParticle * channel80[] ={Neutron,SigmaP,PiPlus1,PiPlus2,PiMinus1,K0S};
    c->AddChannel(width_channel_107,6,channel80);                     // include decay modes
    //PParticle * channel81[] ={Proton1,Proton2,PiPlus1,PiPlus2,PiPlus3,PiPlus4,PiMinus1,PiMinus2,PiMinus3,PiMinus4};
    //c->AddChannel(width_channel_108,10,channel81);                    // include decay modes
    PParticle * channel82[] ={Proton1,Proton2,PiPlus1,PiMinus1,K0S,K0S};
    c->AddChannel(width_channel_109,6,channel82);                     // include decay modes
#endif    


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TH1F histo("histo", "Lambda1405 mass", 100, 1.2, 1.6);
 
    //p_p->Do(&histo,"_x = [Lambda1405]->M()");
    //p_p->Do("echo -------------------------------------------");
    //p_p->Do("foreach(*) ; [*]->Print()");



    p_p->InitReaction(q,c);              // initialize the reaction
    p_p->Print();

    p_p->loop(nEvents,0,"strange",1,0,1,0,0);



}
