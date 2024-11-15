use super::{
    super::super::test::{assert_eq, TestError},
    Tensor, TensorRank0, TensorRank1, TensorRank1List, TensorRank2, TensorRank2List2D
};

fn get_array() -> [[TensorRank0; 3]; 8] {
    [
        [5.0, 0.0, 0.0],
        [5.0, 5.0, 6.0],
        [3.0, 1.0, 4.0],
        [3.0, 4.0, 2.0],
        [1.0, 0.0, 3.0],
        [1.0, 3.0, 1.0],
        [1.0, 6.0, 0.0],
        [1.0, 1.0, 1.0],
    ]
}

fn get_tensor_rank_1_list() -> TensorRank1List<3, 1, 8> {
    TensorRank1List::new(get_array())
}

fn get_other_tensor_rank_1_list() -> TensorRank1List<3, 1, 8> {
    TensorRank1List::new([
        [3.0, 3.0, 6.0],
        [2.0, 4.0, 3.0],
        [6.0, 2.0, 5.0],
        [5.0, 2.0, 5.0],
        [4.0, 7.0, 2.0],
        [2.0, 6.0, 6.0],
        [2.0, 3.0, 2.0],
        [3.0, 7.0, 5.0],
    ])
}

fn get_tensor_rank_1_list_add_other_tensor_rank_1_list() -> TensorRank1List<3, 1, 8> {
    TensorRank1List::new([
        [8.0, 3.0, 6.0],
        [7.0, 9.0, 9.0],
        [9.0, 3.0, 9.0],
        [8.0, 6.0, 7.0],
        [5.0, 7.0, 5.0],
        [3.0, 9.0, 7.0],
        [3.0, 9.0, 2.0],
        [4.0, 8.0, 6.0],
    ])
}

fn get_tensor_rank_1_list_sub_other_tensor_rank_1_list() -> TensorRank1List<3, 1, 8> {
    TensorRank1List::new([
        [2.0, -3.0, -6.0],
        [3.0, 1.0, 3.0],
        [-3.0, -1.0, -1.0],
        [-2.0, 2.0, -3.0],
        [-3.0, -7.0, 1.0],
        [-1.0, -3.0, -5.0],
        [-1.0, 3.0, -2.0],
        [-2.0, -6.0, -4.0],
    ])
}

fn get_tensor_rank_1_list_mul_other_tensor_rank_1_list() -> TensorRank2<3, 1, 1> {
    TensorRank2::new([[69.0, 70.0, 90.0], [57.0, 73.0, 75.0], [63.0, 70.0, 65.0]])
}

fn get_tensor_rank_2_list_2d() -> TensorRank2List2D<3, 1, 0, 8, 8> {
    TensorRank2List2D::new([[[[ 0.1072107436108023   , -0.04342320055090676  ,
        -0.1871804822852251   ],
       [-0.35108142268659215  , -0.07639925253997215  ,
         0.35692195150593986  ],
       [-0.4622046285591005   ,  0.10029376180801064  ,
        -0.39442011202413496  ]],

      [[-0.2539380078555282   , -0.3650142366680812   ,
         0.07437194399749458  ],
       [-0.1319322994915889   ,  0.27168664191770486  ,
        -0.4911990344154312   ],
       [-0.03365557414108944  , -0.14559082917258848  ,
        -0.4245742405743478   ]],

      [[ 0.19099062106787557  ,  0.19173460392506847  ,
         0.3090950026135688   ],
       [-0.16275818219475569  , -0.014953042310487663 ,
        -0.10443116184652934  ],
       [-0.44708770971313627  , -0.4749535599210014   ,
        -0.26068986174617625  ]],

      [[ 0.019557304154016908 , -0.2763247843066866   ,
         0.2775576868522365   ],
       [-0.33519596535484475  ,  0.35892876844129373  ,
        -0.25917990219012943  ],
       [ 0.2882794451585996   ,  0.24198409224028095  ,
        -0.06921464946213107  ]],

      [[-0.17811221285495693  ,  0.4138264199484897   ,
        -0.4764764195400887   ],
       [ 0.42834560917302744  ,  0.09189889688521424  ,
         0.06125584260371364  ],
       [ 0.26826017153100423  , -0.35656150965763334  ,
         0.4414374112639552   ]],

      [[ 0.3422050742861925   ,  0.2131217413230616   ,
        -0.3456801019363118   ],
       [-0.11336695338214775  ,  0.34068785878441044  ,
        -0.3798108913989221   ],
       [-0.10281563052490172  ,  0.20468548426219013  ,
        -0.3241772798727708   ]],

      [[-0.19694659782094437  , -0.08849138242366172  ,
        -0.3340865922495876   ],
       [ 0.07711393832100899  , -0.464692614450238    ,
         0.37417950618853446  ],
       [-0.47683488263525686  , -0.2079155502413318   ,
        -0.4494310239119663   ]],

      [[ 0.052297609442964355 , -0.10927254463218117  ,
         0.14988298543325784  ],
       [-0.1264510302352716   ,  0.0308589933285931   ,
        -0.2370088523375995   ],
       [ 0.443608201301273    , -0.014652018583727111 ,
        -0.49923497988194565  ]]],


     [[[-0.2021591759185496   , -0.3240712400085568   ,
         0.13767559471028168  ],
       [-0.3077563727993886   , -0.13592627772240873  ,
         0.3975041444837816   ],
       [ 0.11398508136122931  ,  0.1420090006052095   ,
        -0.034464704105541144 ]],

      [[-0.010443688227337433 ,  0.15370098498126328  ,
        -0.012815581893079298 ],
       [ 0.48533288608396796  ,  0.20984653275296583  ,
         0.3904989365232152   ],
       [ 0.08704631542160624  , -0.41205094213837257  ,
        -0.22669085695022106  ]],

      [[-0.009522178202827414 , -0.13485131128213712  ,
        -0.3940081108027941   ],
       [-0.2048449971811962   , -0.3546253608542237   ,
        -0.4743171598456104   ],
       [-0.05983130468800768  , -0.2386770248769039   ,
         0.004136367272426167 ]],

      [[-0.2755925442749855   , -0.45094541822813194  ,
         0.3795894677684446   ],
       [ 0.17084038359020282  ,  0.3607005547173684   ,
         0.1479415592333948   ],
       [ 0.31755780246742693  , -0.051817093482198984 ,
        -0.322769628763529    ]],

      [[ 0.182746186273829    ,  0.3081661589332926   ,
         0.3810835140433062   ],
       [ 0.10627987750197865  ,  0.2691325655691511   ,
         0.052165771122818505 ],
       [-0.49474418597820036  ,  0.20150774710962616  ,
        -0.23630990561271525  ]],

      [[ 0.41882519896534376  , -0.04675718255742234  ,
         0.25823972208792445  ],
       [-0.31466187893426034  , -0.004119734186210922 ,
         0.4396711715042375   ],
       [ 0.4669712321466173   ,  0.2751760894485914   ,
         0.10664546713095846  ]],

      [[-0.4872532554053184   ,  0.2122611665750268   ,
        -0.3072950827827039   ],
       [-0.48495157037017655  , -0.4457014527399976   ,
         0.4402177803528132   ],
       [ 0.40213709543474174  ,  0.15361666862804157  ,
         0.2758780847869903   ]],

      [[ 0.27122063477020797  ,  0.2947660051341414   ,
         0.29460156360002854  ],
       [ 0.06401104480886077  , -0.3577261235076473   ,
         0.12116485664082166  ],
       [-0.48866328981374063  , -0.0944367621680483   ,
         0.19329782797707407  ]]],


     [[[ 0.19062722770049112  , -0.08759884436937015  ,
         0.33106260291423806  ],
       [-0.2977978740060243   , -0.32567885264981344  ,
         0.08398508381575676  ],
       [-0.4395161761351305   , -0.3751672991413765   ,
         0.38403802755636995  ]],

      [[-0.37278368943001317  , -0.49533245102184464  ,
        -0.09413735698754322  ],
       [ 0.21721460767368495  ,  0.17759990156493766  ,
         0.19107939602920787  ],
       [ 0.3981516239444868   , -0.3525799891044855   ,
         0.3525079731042643   ]],

      [[-0.4284100888109854   ,  0.07055268304148943  ,
         0.4216975914938028   ],
       [ 0.4795459304603278   ,  0.17394192827386512  ,
        -0.32108481693909674  ],
       [ 0.3467957307244245   , -0.29442538122721196  ,
        -0.2696369312634984   ]],

      [[ 0.025044967450541522 ,  0.10404118347126023  ,
        -0.21812339615536014  ],
       [-0.041472435935500984 ,  0.46880509050341024  ,
         0.17925107319863098  ],
       [ 0.1998370714391715   ,  0.45115831011023366  ,
         0.21436926044134175  ]],

      [[-0.3562104097003427   ,  0.4766299475209881   ,
        -0.08362420821371175  ],
       [ 0.1895448815197316   , -0.08862774811931273  ,
        -0.11266427660791578  ],
       [ 0.1994781647275532   , -0.11072616226268295  ,
         0.4310329284787696   ]],

      [[ 0.1537526719608856   , -0.014854206337902687 ,
         0.04871837946238755  ],
       [ 0.1538614069679196   , -0.02399754807520582  ,
         0.24089614539528437  ],
       [-0.42657693397459173  , -0.2277693843428067   ,
        -0.19367647377849173  ]],

      [[-0.31689751406979905  , -0.4315317399612447   ,
        -0.18018070100992023  ],
       [-0.4485988179220173   ,  0.29601810614319013  ,
        -0.4645272677736719   ],
       [ 0.24747700511273174  ,  0.2868859253388303   ,
        -0.26280619247459414  ]],

      [[ 0.14663345683780993  ,  0.08728180234867078  ,
         0.09756784365610582  ],
       [-0.13175338221388222  , -0.2539780532704241   ,
        -0.13695453915271694  ],
       [-0.4596067103724808   , -0.07622541731280219  ,
         0.10979119955700323  ]]],


     [[[-0.10120171528285138  ,  0.2339842745058779   ,
        -0.26678407794630843  ],
       [ 0.16737415841355663  , -0.18716647781629858  ,
         0.47631093649070655  ],
       [-0.4056069416064043   , -0.23436487084086943  ,
        -0.39161058296315887  ]],

      [[-0.45265570785421616  , -0.4174408784860676   ,
         0.03585070167364535  ],
       [-0.3074549097126109   , -0.43303827875448886  ,
        -0.14344085672539497  ],
       [ 0.3220674681581456   , -0.36622743704983696  ,
        -0.018793205291239734 ]],

      [[-0.338227286761532    ,  0.30593084137910775  ,
        -0.33327563010464234  ],
       [-0.4072820733788014   , -0.04836866202293644  ,
        -0.4135112255363158   ],
       [-0.4466167631337985   , -0.025509103977717218 ,
         0.11243348248769036  ]],

      [[ 0.06941320674382245  , -0.052437683721715866 ,
        -0.358556904940898    ],
       [-0.1953087568080255   ,  0.18231002545772634  ,
        -0.18369349039268423  ],
       [-0.056324954251208315 , -0.3768253237409638   ,
        -0.43377873928609634  ]],

      [[ 0.016022477857021156 , -0.21032236934333504  ,
        -0.2741566756053633   ],
       [-0.2529350462961769   , -0.22727053688003718  ,
         0.16172395901517012  ],
       [-0.2862591088754758   , -0.1941888897571824   ,
        -0.3253288533981188   ]],

      [[-0.20616639955766636  ,  0.47517241431800505  ,
        -0.13794450101503575  ],
       [ 0.3334987251024907   ,  0.007445169675950547 ,
         0.11214180209535296  ],
       [-0.2851886266429289   , -0.33185347960456857  ,
         0.22961075278335608  ]],

      [[ 0.38727892098036476  , -0.23125136277388536  ,
         0.14704033316845133  ],
       [ 0.28244017367799135  , -0.344612822458487    ,
        -0.2579347205175998   ],
       [ 0.12375207784787412  ,  0.4031812341663179   ,
         0.23336470171840817  ]],

      [[ 0.29344382184232853  , -0.49020378793196573  ,
         0.4184238480155871   ],
       [ 0.4264164620510238   , -0.25737624143849835  ,
         0.40277938991398066  ],
       [-0.21664147812318013  , -0.29679858729369146  ,
        -0.21245077320729866  ]]],


     [[[-0.10773683564964309  ,  0.45502891146345326  ,
        -0.39785437995449147  ],
       [ 0.1211926453681399   , -0.15016566105120244  ,
         0.3307509464675549   ],
       [ 0.051523912923697734 ,  0.4081842358078295   ,
        -0.19361808545965398  ]],

      [[-0.330860689915487    ,  0.24950431237265336  ,
        -0.05050208812686041  ],
       [-0.4130733034169801   , -0.27946569034525537  ,
        -0.47362413571094764  ],
       [ 0.19259045231074423  , -0.15566434823405495  ,
         0.014029656717868044 ]],

      [[ 0.24700488467241632  ,  0.15223735511771852  ,
        -0.3616972614773778   ],
       [-0.12155880469623137  , -0.00911705183942646  ,
         0.331421045074234    ],
       [ 0.2816613566423991   , -0.25654182728487207  ,
         0.3558385591090334   ]],

      [[ 0.27390285445671314  ,  0.2346734591322901   ,
        -0.21702912102838812  ],
       [ 0.0271342328784574   , -0.2767102984744013   ,
        -0.45670699788773217  ],
       [-0.4049550213338031   ,  0.2021594790861585   ,
         0.464774317013659    ]],

      [[ 0.23363511432114703  ,  0.004931874426084115 ,
         0.15500642673614373  ],
       [ 0.0046766496227130805, -0.03112678179344186  ,
         0.3479221234006099   ],
       [-0.18651169357330133  , -0.2090453689153975   ,
        -0.1668828568242362   ]],

      [[-0.2453390174661232   , -0.14318796932509115  ,
        -0.15574132191258627  ],
       [ 0.3220986404033934   ,  0.38461696470643136  ,
        -0.4696300235384683   ],
       [ 0.13351185378878871  ,  0.05659751074612884  ,
         0.3005252319840386   ]],

      [[-0.4621096307718717   ,  0.037008750019236425 ,
         0.16926336111926532  ],
       [ 0.34038431705115324  ,  0.00798059195686196  ,
        -0.09079640303839076  ],
       [ 0.4470547552489711   ,  0.1604563600454757   ,
        -0.24428320007259852  ]],

      [[-0.08417178358838429  ,  0.39315967754185666  ,
         0.25406896764499387  ],
       [ 0.31381383366661963  ,  0.2843544947654961   ,
        -0.3688175039763809   ],
       [-0.2255479745030894   ,  0.17507407791503793  ,
        -0.3028062402801155   ]]],


     [[[ 0.015717412247830964 ,  0.22720039175809048  ,
        -0.09242702668852576  ],
       [ 0.3595972787376489   ,  0.2686525966108094   ,
         0.24239367731171668  ],
       [ 0.12315768775498503  ,  0.3652090941256053   ,
        -0.46706197083529055  ]],

      [[ 0.18664746833385215  , -0.06700500013776589  ,
        -0.37744778578532034  ],
       [-0.19232434472700943  ,  0.048708811822774645 ,
        -0.2344497205548972   ],
       [ 0.3758706340944552   , -0.24199109293513288  ,
         0.02959693852786438  ]],

      [[ 0.10313258168159445  , -0.19550861582085766  ,
         0.35988620709461105  ],
       [ 0.1955917709740873   ,  0.3161459441103066   ,
        -0.4535761720300354   ],
       [-0.3789788154744902   ,  0.09329346043266296  ,
        -0.3563522018055668   ]],

      [[ 0.3218795459162863   , -0.22087926766709232  ,
        -0.45407690331995854  ],
       [ 0.19231357212455535  ,  0.28891307797369103  ,
        -0.37925458853974103  ],
       [ 0.11016169558037614  ,  0.056080088274753415 ,
        -0.4843846447621023   ]],

      [[-0.23075485292560283  , -0.39827755644433627  ,
         0.26024644617575865  ],
       [-0.3773777193824285   , -0.13348875213833833  ,
        -0.4315926406851265   ],
       [ 0.03734162418659037  ,  0.19666306848271065  ,
        -0.3012303202280474   ]],

      [[ 0.47433597241457315  ,  0.4568268025661991   ,
        -0.29597226549932687  ],
       [-0.03708338506689046  ,  0.11760508737043684  ,
         0.04654823656267748  ],
       [ 0.27159598942318575  , -0.30513529938139916  ,
         0.11651724505306726  ]],

      [[-0.49284830577580285  , -0.30163237022836387  ,
         0.10366631803541904  ],
       [-0.010632491900066432 , -0.4080255270715083   ,
         0.35617296132649967  ],
       [-0.09462715917895193  ,  0.03612290071310653  ,
        -0.22275932897764084  ]],

      [[ 0.17316471160729485  ,  0.07880456455653817  ,
        -0.3000505755299464   ],
       [ 0.31458292775059227  , -0.13234543511201358  ,
         0.14489370549033198  ],
       [-0.40775534889495013  ,  0.46336330918295954  ,
        -0.4462789288195016   ]]],


     [[[ 0.09007265486265892  ,  0.30427318757657273  ,
        -0.24451724136753317  ],
       [ 0.4036479587909003   , -0.39859741747252564  ,
         0.42744038825774266  ],
       [ 0.18640446038621106  ,  0.08727657680251955  ,
         0.2551390368624138   ]],

      [[ 0.1419871583656348   ,  0.19627216854582008  ,
        -0.15875257132898712  ],
       [ 0.44308807028553676  , -0.3772553944262872   ,
         0.11168320454398206  ],
       [ 0.488605367437576    ,  0.4425631371535901   ,
        -0.1111480248921548   ]],

      [[ 0.3503659779962107   , -0.2303269082563787   ,
         0.44296508988413297  ],
       [ 0.23050562081696668  ,  0.19500799997712348  ,
        -0.07356013976360198  ],
       [ 0.26464932708670286  , -0.4714366673535958   ,
         0.37401764333953913  ]],

      [[-0.45523469289239005  ,  0.32363042233777284  ,
        -0.3597612923806002   ],
       [ 0.4176692744955405   ,  0.07523062778735001  ,
        -0.34808369725596056  ],
       [ 0.27540252194070525  , -0.3814219048747137   ,
         0.220534961130619    ]],

      [[ 0.25281219441990066  ,  0.36409673604911075  ,
         0.4177484003247197   ],
       [ 0.1254388647216811   ,  0.4216966791925947   ,
        -0.060220429107027074 ],
       [ 0.19806734902117662  ,  0.14813491885217012  ,
        -0.31919298227026327  ]],

      [[-0.2539909533508231   , -0.010678372219592447 ,
        -0.3215853047475623   ],
       [-0.4624454628599439   , -0.2635556052597601   ,
        -0.3514984738394491   ],
       [ 0.33687393148960565  , -0.4040504327456197   ,
         0.47468809194969197  ]],

      [[ 0.4809886829350227   ,  0.4697563109136348   ,
         0.20029270173370728  ],
       [-0.434549599600643    , -0.14296948545290133  ,
        -0.33098599292515296  ],
       [-0.27096903860623334  , -0.2528267422379076   ,
        -0.475464713235747    ]],

      [[-0.47750865175014834  ,  0.04580328542372292  ,
         0.41845951713939267  ],
       [-0.285564306662332    , -0.030806026399835207 ,
        -0.03374793050610003  ],
       [ 0.031158504278914978 ,  0.14468400727305264  ,
        -0.11585441307449551  ]]],


     [[[ 0.32182310154029536  , -0.2745473190723162   ,
         0.48639945920653194  ],
       [-0.4793631911651609   , -0.3333548973436097   ,
         0.13461353355252725  ],
       [-0.4842840970840331   , -0.3323889455835023   ,
         0.29733951821844595  ]],

      [[ 0.17220562799318206  , -0.3290229336905529   ,
         0.3045034624345301   ],
       [-0.4690516699274022   , -0.195901507119286    ,
         0.3197523058058891   ],
       [ 0.130023594031362    ,  0.009877450530229814 ,
         0.04566186102919256  ]],

      [[ 0.2226787722990089   ,  0.06363266453927385  ,
         0.3373532580869739   ],
       [-0.21569876477718408  , -0.008816713637215479 ,
         0.20035988600827714  ],
       [ 0.22350754317814747  ,  0.04818665144190648  ,
         0.18192190169388978  ]],

      [[ 0.48063988499941346  ,  0.14106224001336753  ,
        -0.24080768563120147  ],
       [ 0.26647632188178505  , -0.07522364444074037  ,
         0.32681006911058386  ],
       [ 0.33835996171688654  , -0.2554845133427329   ,
         0.42293945359140417  ]],

      [[ 0.06365268441252425  , -0.352702199416716    ,
        -0.473919023735638    ],
       [ 0.2409958842646336   ,  0.45016687025490376  ,
         0.43224908615835866  ],
       [-0.2357561838048864   , -0.24971696384583486  ,
        -0.08642472830921932  ]],

      [[ 0.059237814966723645 , -0.20254807338745295  ,
        -0.3180478588033352   ],
       [ 0.013830239337064287 ,  0.4296092319016629   ,
         0.4138171910785582   ],
       [ 0.19900852468571273  ,  0.22444499139612684  ,
        -0.25210040288262625  ]],

      [[-0.3857104179316506   , -0.4371746902461279   ,
        -0.12866082437045334  ],
       [ 0.3118514389493242   ,  0.22704669517608622  ,
        -0.484517214076141    ],
       [-0.4086774635676914   ,  0.3719364936377215   ,
         0.0673117848145729   ]],

      [[-0.40970667977900066  , -0.00579943023197349  ,
         0.4282433037697012   ],
       [-0.1976153000375621   ,  0.009232935205557347 ,
         0.15212194604286067  ],
       [-0.36613942522171017  , -0.2409490741396667   ,
         0.20674716322774622  ]]]])
}

fn get_tensor_rank_1_list_div_tensor_rank_2_list_2d() -> TensorRank1List<3, 0, 8> {
    TensorRank1List::new([
        [2.0464212135154436, 6.827698367888031, 11.056836790412634],
        [-22.947026677779508, 10.574323884914493, 21.750457350457467],
        [8.044208887151221, -35.07133502812548, -18.61369225864542],
        [-20.684268614523628, -19.39220610450351, -19.57545635615623],
        [-23.727347518884883, 1.2066552121452456, -7.7838397852853864],
        [-16.43920097416692, 13.970188500416638, -6.588123940005775],
        [-10.813302431346484, -2.066475959739897, -19.36077096462388],
        [-20.146900991605264, -4.004412133345104, -8.310553065057563]
    ])
}

#[test]
fn add_tensor_rank_1_list_to_self() {
    (get_tensor_rank_1_list() + get_other_tensor_rank_1_list())
        .iter()
        .zip(get_tensor_rank_1_list_add_other_tensor_rank_1_list().iter())
        .for_each(|(tensor_rank_1_list_entry, add_tensor_rank_1_list_entry)| {
            tensor_rank_1_list_entry
                .iter()
                .zip(add_tensor_rank_1_list_entry.iter())
                .for_each(
                    |(tensor_rank_1_list_entry_i, add_tensor_rank_1_list_entry_i)| {
                        assert_eq!(tensor_rank_1_list_entry_i, add_tensor_rank_1_list_entry_i)
                    },
                )
        });
}

#[test]
fn add_tensor_rank_1_list_ref_to_self() {
    (get_tensor_rank_1_list() + &get_other_tensor_rank_1_list())
        .iter()
        .zip(get_tensor_rank_1_list_add_other_tensor_rank_1_list().iter())
        .for_each(|(tensor_rank_1_list_entry, add_tensor_rank_1_list_entry)| {
            tensor_rank_1_list_entry
                .iter()
                .zip(add_tensor_rank_1_list_entry.iter())
                .for_each(
                    |(tensor_rank_1_list_entry_i, add_tensor_rank_1_list_entry_i)| {
                        assert_eq!(tensor_rank_1_list_entry_i, add_tensor_rank_1_list_entry_i)
                    },
                )
        });
}

#[test]
fn add_assign_tensor_rank_1_list() {
    let mut tensor_rank_1_list = get_tensor_rank_1_list();
    tensor_rank_1_list += get_other_tensor_rank_1_list();
    tensor_rank_1_list
        .iter()
        .zip(get_tensor_rank_1_list_add_other_tensor_rank_1_list().iter())
        .for_each(|(tensor_rank_1_list_entry, add_tensor_rank_1_list_entry)| {
            tensor_rank_1_list_entry
                .iter()
                .zip(add_tensor_rank_1_list_entry.iter())
                .for_each(
                    |(tensor_rank_1_list_entry_i, add_tensor_rank_1_list_entry_i)| {
                        assert_eq!(tensor_rank_1_list_entry_i, add_tensor_rank_1_list_entry_i)
                    },
                )
        });
}

#[test]
fn add_assign_tensor_rank_1_list_ref() {
    let mut tensor_rank_1_list = get_tensor_rank_1_list();
    tensor_rank_1_list += &get_other_tensor_rank_1_list();
    tensor_rank_1_list
        .iter()
        .zip(get_tensor_rank_1_list_add_other_tensor_rank_1_list().iter())
        .for_each(|(tensor_rank_1_list_entry, add_tensor_rank_1_list_entry)| {
            tensor_rank_1_list_entry
                .iter()
                .zip(add_tensor_rank_1_list_entry.iter())
                .for_each(
                    |(tensor_rank_1_list_entry_i, add_tensor_rank_1_list_entry_i)| {
                        assert_eq!(tensor_rank_1_list_entry_i, add_tensor_rank_1_list_entry_i)
                    },
                )
        });
}

#[test]
fn as_array() {
    get_tensor_rank_1_list()
        .as_array()
        .iter()
        .zip(get_array().iter())
        .for_each(|(tensor_rank_1_as_array_entry, array_entry)| {
            tensor_rank_1_as_array_entry
                .iter()
                .zip(array_entry.iter())
                .for_each(|(tensor_rank_1_as_array_entry_i, array_entry_i)| {
                    assert_eq!(tensor_rank_1_as_array_entry_i, array_entry_i)
                })
        });
}

#[test]
fn div_tensor_rank_0_to_self() {
    (get_tensor_rank_1_list() / 3.3)
        .iter()
        .zip(get_array().iter())
        .for_each(|(tensor_rank_1_list_entry, array_entry)| {
            tensor_rank_1_list_entry
                .iter()
                .zip(array_entry.iter())
                .for_each(|(tensor_rank_1_list_entry_i, array_entry_i)| {
                    assert_eq!(tensor_rank_1_list_entry_i, &(array_entry_i / 3.3))
                })
        });
}

#[test]
#[allow(clippy::op_ref)]
fn div_tensor_rank_0_ref_to_self() {
    (get_tensor_rank_1_list() / &3.3)
        .iter()
        .zip(get_array().iter())
        .for_each(|(tensor_rank_1_list_entry, array_entry)| {
            tensor_rank_1_list_entry
                .iter()
                .zip(array_entry.iter())
                .for_each(|(tensor_rank_1_list_entry_i, array_entry_i)| {
                    assert_eq!(tensor_rank_1_list_entry_i, &(array_entry_i / 3.3))
                })
        });
}

#[test]
fn div_assign_tensor_rank_0() {
    let mut tensor_rank_1_list = get_tensor_rank_1_list();
    tensor_rank_1_list /= 3.3;
    tensor_rank_1_list.iter().zip(get_array().iter()).for_each(
        |(tensor_rank_1_list_entry, array_entry)| {
            tensor_rank_1_list_entry
                .iter()
                .zip(array_entry.iter())
                .for_each(|(tensor_rank_1_list_entry_i, array_entry_i)| {
                    assert_eq!(tensor_rank_1_list_entry_i, &(array_entry_i / 3.3))
                })
        },
    );
}

#[test]
fn div_assign_tensor_rank_0_ref() {
    let mut tensor_rank_1_list = get_tensor_rank_1_list();
    tensor_rank_1_list /= &3.3;
    tensor_rank_1_list.iter().zip(get_array().iter()).for_each(
        |(tensor_rank_1_list_entry, array_entry)| {
            tensor_rank_1_list_entry
                .iter()
                .zip(array_entry.iter())
                .for_each(|(tensor_rank_1_list_entry_i, array_entry_i)| {
                    assert_eq!(tensor_rank_1_list_entry_i, &(array_entry_i / 3.3))
                })
        },
    );
}

#[test]
fn div_tensor_rank_4_to_self_dim_3() -> Result<(), TestError> {
    assert_eq(
        &(get_tensor_rank_1_list() / get_tensor_rank_2_list_2d()),
        &get_tensor_rank_1_list_div_tensor_rank_2_list_2d(),
    )
}

#[test]
fn from_iter() {
    let into_iterator = get_tensor_rank_1_list().0;
    let tensor_rank_1_list = TensorRank1List::<3, 1, 8>::from_iter(get_tensor_rank_1_list().0);
    tensor_rank_1_list
        .iter()
        .zip(into_iterator)
        .for_each(|(tensor_rank_1_list_entry, entry)| {
            tensor_rank_1_list_entry.iter().zip(entry.iter()).for_each(
                |(tensor_rank_1_list_entry_i, entry_i)| {
                    assert_eq!(tensor_rank_1_list_entry_i, entry_i)
                },
            )
        });
}

#[test]
fn iter() {
    get_tensor_rank_1_list()
        .iter()
        .zip(get_array().iter())
        .for_each(|(tensor_rank_1_entry, array_entry)| {
            tensor_rank_1_entry.iter().zip(array_entry.iter()).for_each(
                |(tensor_rank_1_entry_i, array_entry_i)| {
                    assert_eq!(tensor_rank_1_entry_i, array_entry_i)
                },
            )
        });
}

#[test]
fn iter_mut() {
    get_tensor_rank_1_list()
        .iter_mut()
        .zip(get_array().iter_mut())
        .for_each(|(tensor_rank_1_entry, array_entry)| {
            tensor_rank_1_entry
                .iter_mut()
                .zip(array_entry.iter_mut())
                .for_each(|(tensor_rank_1_entry_i, array_entry_i)| {
                    assert_eq!(tensor_rank_1_entry_i, array_entry_i)
                })
        });
}

#[test]
fn mul_tensor_rank_0_to_self() {
    (get_tensor_rank_1_list() * 3.3)
        .iter()
        .zip(get_array().iter())
        .for_each(|(tensor_rank_1_list_entry, array_entry)| {
            tensor_rank_1_list_entry
                .iter()
                .zip(array_entry.iter())
                .for_each(|(tensor_rank_1_list_entry_i, array_entry_i)| {
                    assert_eq!(tensor_rank_1_list_entry_i, &(array_entry_i * 3.3))
                })
        });
}

#[test]
#[allow(clippy::op_ref)]
fn mul_tensor_rank_0_ref_to_self() {
    (get_tensor_rank_1_list() * &3.3)
        .iter()
        .zip(get_array().iter())
        .for_each(|(tensor_rank_1_list_entry, array_entry)| {
            tensor_rank_1_list_entry
                .iter()
                .zip(array_entry.iter())
                .for_each(|(tensor_rank_1_list_entry_i, array_entry_i)| {
                    assert_eq!(tensor_rank_1_list_entry_i, &(array_entry_i * 3.3))
                })
        });
}

#[test]
#[allow(clippy::op_ref)]
fn mul_tensor_rank_0_ref_to_self_ref() {
    (&get_tensor_rank_1_list() * &3.3)
        .iter()
        .zip(get_array().iter())
        .for_each(|(tensor_rank_1_list_entry, array_entry)| {
            tensor_rank_1_list_entry
                .iter()
                .zip(array_entry.iter())
                .for_each(|(tensor_rank_1_list_entry_i, array_entry_i)| {
                    assert_eq!(tensor_rank_1_list_entry_i, &(array_entry_i * 3.3))
                })
        });
}

#[test]
fn mul_assign_tensor_rank_0() {
    let mut tensor_rank_1_list = get_tensor_rank_1_list();
    tensor_rank_1_list *= 3.3;
    tensor_rank_1_list.iter().zip(get_array().iter()).for_each(
        |(tensor_rank_1_list_entry, array_entry)| {
            tensor_rank_1_list_entry
                .iter()
                .zip(array_entry.iter())
                .for_each(|(tensor_rank_1_list_entry_i, array_entry_i)| {
                    assert_eq!(tensor_rank_1_list_entry_i, &(array_entry_i * 3.3))
                })
        },
    );
}

#[test]
fn mul_assign_tensor_rank_0_ref() {
    let mut tensor_rank_1_list = get_tensor_rank_1_list();
    tensor_rank_1_list *= &3.3;
    tensor_rank_1_list.iter().zip(get_array().iter()).for_each(
        |(tensor_rank_1_list_entry, array_entry)| {
            tensor_rank_1_list_entry
                .iter()
                .zip(array_entry.iter())
                .for_each(|(tensor_rank_1_list_entry_i, array_entry_i)| {
                    assert_eq!(tensor_rank_1_list_entry_i, &(array_entry_i * 3.3))
                })
        },
    );
}

#[test]
fn mul_tensor_rank_1_list_to_self() {
    (get_tensor_rank_1_list() * get_other_tensor_rank_1_list())
        .iter()
        .zip(get_tensor_rank_1_list_mul_other_tensor_rank_1_list().iter())
        .for_each(|(tensor_rank_1_list_entry, mul_entry)| {
            tensor_rank_1_list_entry
                .iter()
                .zip(mul_entry.iter())
                .for_each(|(tensor_rank_1_list_entry_i, mul_entry_i)| {
                    assert_eq!(tensor_rank_1_list_entry_i, mul_entry_i)
                })
        });
}

#[test]
fn mul_tensor_rank_1_list_ref_to_self() {
    (get_tensor_rank_1_list() * &get_other_tensor_rank_1_list())
        .iter()
        .zip(get_tensor_rank_1_list_mul_other_tensor_rank_1_list().iter())
        .for_each(|(tensor_rank_1_list_entry, mul_entry)| {
            tensor_rank_1_list_entry
                .iter()
                .zip(mul_entry.iter())
                .for_each(|(tensor_rank_1_list_entry_i, mul_entry_i)| {
                    assert_eq!(tensor_rank_1_list_entry_i, mul_entry_i)
                })
        });
}

#[test]
fn mul_tensor_rank_1_list_to_self_ref() {
    (&get_tensor_rank_1_list() * get_other_tensor_rank_1_list())
        .iter()
        .zip(get_tensor_rank_1_list_mul_other_tensor_rank_1_list().iter())
        .for_each(|(tensor_rank_1_list_entry, mul_entry)| {
            tensor_rank_1_list_entry
                .iter()
                .zip(mul_entry.iter())
                .for_each(|(tensor_rank_1_list_entry_i, mul_entry_i)| {
                    assert_eq!(tensor_rank_1_list_entry_i, mul_entry_i)
                })
        });
}

#[test]
fn mul_tensor_rank_1_list_ref_to_self_ref() {
    (&get_tensor_rank_1_list() * &get_other_tensor_rank_1_list())
        .iter()
        .zip(get_tensor_rank_1_list_mul_other_tensor_rank_1_list().iter())
        .for_each(|(tensor_rank_1_list_entry, mul_entry)| {
            tensor_rank_1_list_entry
                .iter()
                .zip(mul_entry.iter())
                .for_each(|(tensor_rank_1_list_entry_i, mul_entry_i)| {
                    assert_eq!(tensor_rank_1_list_entry_i, mul_entry_i)
                })
        });
}

#[test]
fn new() {
    get_tensor_rank_1_list()
        .iter()
        .zip(get_array().iter())
        .for_each(|(tensor_rank_1_entry, array_entry)| {
            tensor_rank_1_entry.iter().zip(array_entry.iter()).for_each(
                |(tensor_rank_1_entry_i, array_entry_i)| {
                    assert_eq!(tensor_rank_1_entry_i, array_entry_i)
                },
            )
        });
}

#[test]
fn size() {
    assert_eq!(
        std::mem::size_of::<TensorRank1List::<3, 1, 8>>(),
        std::mem::size_of::<[TensorRank1::<3, 1>; 8]>()
    )
}

#[test]
fn sub_tensor_rank_1_list_to_self() {
    (get_tensor_rank_1_list() - get_other_tensor_rank_1_list())
        .iter()
        .zip(get_tensor_rank_1_list_sub_other_tensor_rank_1_list().iter())
        .for_each(|(tensor_rank_1_list_entry, add_tensor_rank_1_list_entry)| {
            tensor_rank_1_list_entry
                .iter()
                .zip(add_tensor_rank_1_list_entry.iter())
                .for_each(
                    |(tensor_rank_1_list_entry_i, add_tensor_rank_1_list_entry_i)| {
                        assert_eq!(tensor_rank_1_list_entry_i, add_tensor_rank_1_list_entry_i)
                    },
                )
        });
}

#[test]
fn sub_tensor_rank_1_list_ref_to_self() {
    (get_tensor_rank_1_list() - &get_other_tensor_rank_1_list())
        .iter()
        .zip(get_tensor_rank_1_list_sub_other_tensor_rank_1_list().iter())
        .for_each(|(tensor_rank_1_list_entry, add_tensor_rank_1_list_entry)| {
            tensor_rank_1_list_entry
                .iter()
                .zip(add_tensor_rank_1_list_entry.iter())
                .for_each(
                    |(tensor_rank_1_list_entry_i, add_tensor_rank_1_list_entry_i)| {
                        assert_eq!(tensor_rank_1_list_entry_i, add_tensor_rank_1_list_entry_i)
                    },
                )
        });
}

#[test]
fn sub_assign_tensor_rank_1_list() {
    let mut tensor_rank_1_list = get_tensor_rank_1_list();
    tensor_rank_1_list -= get_other_tensor_rank_1_list();
    tensor_rank_1_list
        .iter()
        .zip(get_tensor_rank_1_list_sub_other_tensor_rank_1_list().iter())
        .for_each(|(tensor_rank_1_list_entry, add_tensor_rank_1_list_entry)| {
            tensor_rank_1_list_entry
                .iter()
                .zip(add_tensor_rank_1_list_entry.iter())
                .for_each(
                    |(tensor_rank_1_list_entry_i, add_tensor_rank_1_list_entry_i)| {
                        assert_eq!(tensor_rank_1_list_entry_i, add_tensor_rank_1_list_entry_i)
                    },
                )
        });
}

#[test]
fn sub_assign_tensor_rank_1_list_ref() {
    let mut tensor_rank_1_list = get_tensor_rank_1_list();
    tensor_rank_1_list -= &get_other_tensor_rank_1_list();
    tensor_rank_1_list
        .iter()
        .zip(get_tensor_rank_1_list_sub_other_tensor_rank_1_list().iter())
        .for_each(|(tensor_rank_1_list_entry, add_tensor_rank_1_list_entry)| {
            tensor_rank_1_list_entry
                .iter()
                .zip(add_tensor_rank_1_list_entry.iter())
                .for_each(
                    |(tensor_rank_1_list_entry_i, add_tensor_rank_1_list_entry_i)| {
                        assert_eq!(tensor_rank_1_list_entry_i, add_tensor_rank_1_list_entry_i)
                    },
                )
        });
}

#[test]
fn zero() {
    TensorRank1List::<3, 1, 8>::zero()
        .iter()
        .for_each(|tensor_rank_1_entry| {
            tensor_rank_1_entry
                .iter()
                .for_each(|tensor_rank_1_entry_i| assert_eq!(tensor_rank_1_entry_i, &0.0))
        });
}
