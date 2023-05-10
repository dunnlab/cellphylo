#!/bin/bash

#run_scjack_parallel_540_TBE.sh

#run Jumble in parallel
#only want to set up directories once
#need 1 empty outfile and outtree so you can indicate indiividual names for output files

cd /contml

mkdir outtrees_scjack_540_TBE
mkdir outfiles_scjack_540_TBE/

#make an empty outfile and outtree so contml prompts for custom file names
touch outfile
touch outtree

#ensure that there is no 'infile'. Will throw up error if already no infile.
#rm infile

echo "scjack jumble search preliminaries complete"

#send in the command. Remove last `&`. #chatgpt

#for i in {1..50}; do printf "./script_scjack_540_TBE_$i.sh & %.0s"; done

#run the script 50 times in parallel

./script_scjack_540_TBE_1.sh & ./script_scjack_540_TBE_2.sh & ./script_scjack_540_TBE_3.sh & ./script_scjack_540_TBE_4.sh & ./script_scjack_540_TBE_5.sh & ./script_scjack_540_TBE_6.sh & ./script_scjack_540_TBE_7.sh & ./script_scjack_540_TBE_8.sh & ./script_scjack_540_TBE_9.sh & ./script_scjack_540_TBE_10.sh & ./script_scjack_540_TBE_11.sh & ./script_scjack_540_TBE_12.sh & ./script_scjack_540_TBE_13.sh & ./script_scjack_540_TBE_14.sh & ./script_scjack_540_TBE_15.sh & ./script_scjack_540_TBE_16.sh & ./script_scjack_540_TBE_17.sh & ./script_scjack_540_TBE_18.sh & ./script_scjack_540_TBE_19.sh & ./script_scjack_540_TBE_20.sh & ./script_scjack_540_TBE_21.sh & ./script_scjack_540_TBE_22.sh & ./script_scjack_540_TBE_23.sh & ./script_scjack_540_TBE_24.sh & ./script_scjack_540_TBE_25.sh & ./script_scjack_540_TBE_26.sh & ./script_scjack_540_TBE_27.sh & ./script_scjack_540_TBE_28.sh & ./script_scjack_540_TBE_29.sh & ./script_scjack_540_TBE_30.sh & ./script_scjack_540_TBE_31.sh & ./script_scjack_540_TBE_32.sh & ./script_scjack_540_TBE_33.sh & ./script_scjack_540_TBE_34.sh & ./script_scjack_540_TBE_35.sh & ./script_scjack_540_TBE_36.sh & ./script_scjack_540_TBE_37.sh & ./script_scjack_540_TBE_38.sh & ./script_scjack_540_TBE_39.sh & ./script_scjack_540_TBE_40.sh & ./script_scjack_540_TBE_41.sh & ./script_scjack_540_TBE_42.sh & ./script_scjack_540_TBE_43.sh & ./script_scjack_540_TBE_44.sh & ./script_scjack_540_TBE_45.sh & ./script_scjack_540_TBE_46.sh & ./script_scjack_540_TBE_47.sh & ./script_scjack_540_TBE_48.sh & ./script_scjack_540_TBE_49.sh & ./script_scjack_540_TBE_50.sh

echo "Commands 1-50 complete"

./script_scjack_540_TBE_51.sh & ./script_scjack_540_TBE_52.sh & ./script_scjack_540_TBE_53.sh & ./script_scjack_540_TBE_54.sh & ./script_scjack_540_TBE_55.sh & ./script_scjack_540_TBE_56.sh & ./script_scjack_540_TBE_57.sh & ./script_scjack_540_TBE_58.sh & ./script_scjack_540_TBE_59.sh & ./script_scjack_540_TBE_60.sh & ./script_scjack_540_TBE_61.sh & ./script_scjack_540_TBE_62.sh & ./script_scjack_540_TBE_63.sh & ./script_scjack_540_TBE_64.sh & ./script_scjack_540_TBE_65.sh & ./script_scjack_540_TBE_66.sh & ./script_scjack_540_TBE_67.sh & ./script_scjack_540_TBE_68.sh & ./script_scjack_540_TBE_69.sh & ./script_scjack_540_TBE_70.sh & ./script_scjack_540_TBE_71.sh & ./script_scjack_540_TBE_72.sh & ./script_scjack_540_TBE_73.sh & ./script_scjack_540_TBE_74.sh & ./script_scjack_540_TBE_75.sh & ./script_scjack_540_TBE_76.sh & ./script_scjack_540_TBE_77.sh & ./script_scjack_540_TBE_78.sh & ./script_scjack_540_TBE_79.sh & ./script_scjack_540_TBE_80.sh & ./script_scjack_540_TBE_81.sh & ./script_scjack_540_TBE_82.sh & ./script_scjack_540_TBE_83.sh & ./script_scjack_540_TBE_84.sh & ./script_scjack_540_TBE_85.sh & ./script_scjack_540_TBE_86.sh & ./script_scjack_540_TBE_87.sh & ./script_scjack_540_TBE_88.sh & ./script_scjack_540_TBE_89.sh & ./script_scjack_540_TBE_90.sh & ./script_scjack_540_TBE_91.sh & ./script_scjack_540_TBE_92.sh & ./script_scjack_540_TBE_93.sh & ./script_scjack_540_TBE_94.sh & ./script_scjack_540_TBE_95.sh & ./script_scjack_540_TBE_96.sh & ./script_scjack_540_TBE_97.sh & ./script_scjack_540_TBE_98.sh & ./script_scjack_540_TBE_99.sh & ./script_scjack_540_TBE_100.sh

echo "Commands 51-100 complete"

./script_scjack_540_TBE_101.sh & ./script_scjack_540_TBE_102.sh & ./script_scjack_540_TBE_103.sh & ./script_scjack_540_TBE_104.sh & ./script_scjack_540_TBE_105.sh & ./script_scjack_540_TBE_106.sh & ./script_scjack_540_TBE_107.sh & ./script_scjack_540_TBE_108.sh & ./script_scjack_540_TBE_109.sh & ./script_scjack_540_TBE_110.sh & ./script_scjack_540_TBE_111.sh & ./script_scjack_540_TBE_112.sh & ./script_scjack_540_TBE_113.sh & ./script_scjack_540_TBE_114.sh & ./script_scjack_540_TBE_115.sh & ./script_scjack_540_TBE_116.sh & ./script_scjack_540_TBE_117.sh & ./script_scjack_540_TBE_118.sh & ./script_scjack_540_TBE_119.sh & ./script_scjack_540_TBE_120.sh & ./script_scjack_540_TBE_121.sh & ./script_scjack_540_TBE_122.sh & ./script_scjack_540_TBE_123.sh & ./script_scjack_540_TBE_124.sh & ./script_scjack_540_TBE_125.sh & ./script_scjack_540_TBE_126.sh & ./script_scjack_540_TBE_127.sh & ./script_scjack_540_TBE_128.sh & ./script_scjack_540_TBE_129.sh & ./script_scjack_540_TBE_130.sh & ./script_scjack_540_TBE_131.sh & ./script_scjack_540_TBE_132.sh & ./script_scjack_540_TBE_133.sh & ./script_scjack_540_TBE_134.sh & ./script_scjack_540_TBE_135.sh & ./script_scjack_540_TBE_136.sh & ./script_scjack_540_TBE_137.sh & ./script_scjack_540_TBE_138.sh & ./script_scjack_540_TBE_139.sh & ./script_scjack_540_TBE_140.sh & ./script_scjack_540_TBE_141.sh & ./script_scjack_540_TBE_142.sh & ./script_scjack_540_TBE_143.sh & ./script_scjack_540_TBE_144.sh & ./script_scjack_540_TBE_145.sh & ./script_scjack_540_TBE_146.sh & ./script_scjack_540_TBE_147.sh & ./script_scjack_540_TBE_148.sh & ./script_scjack_540_TBE_149.sh & ./script_scjack_540_TBE_150.sh

echo "Commands 101-150 complete"

./script_scjack_540_TBE_151.sh & ./script_scjack_540_TBE_152.sh & ./script_scjack_540_TBE_153.sh & ./script_scjack_540_TBE_154.sh & ./script_scjack_540_TBE_155.sh & ./script_scjack_540_TBE_156.sh & ./script_scjack_540_TBE_157.sh & ./script_scjack_540_TBE_158.sh & ./script_scjack_540_TBE_159.sh & ./script_scjack_540_TBE_160.sh & ./script_scjack_540_TBE_161.sh & ./script_scjack_540_TBE_162.sh & ./script_scjack_540_TBE_163.sh & ./script_scjack_540_TBE_164.sh & ./script_scjack_540_TBE_165.sh & ./script_scjack_540_TBE_166.sh & ./script_scjack_540_TBE_167.sh & ./script_scjack_540_TBE_168.sh & ./script_scjack_540_TBE_169.sh & ./script_scjack_540_TBE_170.sh & ./script_scjack_540_TBE_171.sh & ./script_scjack_540_TBE_172.sh & ./script_scjack_540_TBE_173.sh & ./script_scjack_540_TBE_174.sh & ./script_scjack_540_TBE_175.sh & ./script_scjack_540_TBE_176.sh & ./script_scjack_540_TBE_177.sh & ./script_scjack_540_TBE_178.sh & ./script_scjack_540_TBE_179.sh & ./script_scjack_540_TBE_180.sh & ./script_scjack_540_TBE_181.sh & ./script_scjack_540_TBE_182.sh & ./script_scjack_540_TBE_183.sh & ./script_scjack_540_TBE_184.sh & ./script_scjack_540_TBE_185.sh & ./script_scjack_540_TBE_186.sh & ./script_scjack_540_TBE_187.sh & ./script_scjack_540_TBE_188.sh & ./script_scjack_540_TBE_189.sh & ./script_scjack_540_TBE_190.sh & ./script_scjack_540_TBE_191.sh & ./script_scjack_540_TBE_192.sh & ./script_scjack_540_TBE_193.sh & ./script_scjack_540_TBE_194.sh & ./script_scjack_540_TBE_195.sh & ./script_scjack_540_TBE_196.sh & ./script_scjack_540_TBE_197.sh & ./script_scjack_540_TBE_198.sh & ./script_scjack_540_TBE_199.sh & ./script_scjack_540_TBE_200.sh 

echo "Commands 151-200 complete"

./script_scjack_540_TBE_201.sh & ./script_scjack_540_TBE_202.sh & ./script_scjack_540_TBE_203.sh & ./script_scjack_540_TBE_204.sh & ./script_scjack_540_TBE_205.sh & ./script_scjack_540_TBE_206.sh & ./script_scjack_540_TBE_207.sh & ./script_scjack_540_TBE_208.sh & ./script_scjack_540_TBE_209.sh & ./script_scjack_540_TBE_210.sh & ./script_scjack_540_TBE_211.sh & ./script_scjack_540_TBE_212.sh & ./script_scjack_540_TBE_213.sh & ./script_scjack_540_TBE_214.sh & ./script_scjack_540_TBE_215.sh & ./script_scjack_540_TBE_216.sh & ./script_scjack_540_TBE_217.sh & ./script_scjack_540_TBE_218.sh & ./script_scjack_540_TBE_219.sh & ./script_scjack_540_TBE_220.sh & ./script_scjack_540_TBE_221.sh & ./script_scjack_540_TBE_222.sh & ./script_scjack_540_TBE_223.sh & ./script_scjack_540_TBE_224.sh & ./script_scjack_540_TBE_225.sh & ./script_scjack_540_TBE_226.sh & ./script_scjack_540_TBE_227.sh & ./script_scjack_540_TBE_228.sh & ./script_scjack_540_TBE_229.sh & ./script_scjack_540_TBE_230.sh & ./script_scjack_540_TBE_231.sh & ./script_scjack_540_TBE_232.sh & ./script_scjack_540_TBE_233.sh & ./script_scjack_540_TBE_234.sh & ./script_scjack_540_TBE_235.sh & ./script_scjack_540_TBE_236.sh & ./script_scjack_540_TBE_237.sh & ./script_scjack_540_TBE_238.sh & ./script_scjack_540_TBE_239.sh & ./script_scjack_540_TBE_240.sh & ./script_scjack_540_TBE_241.sh & ./script_scjack_540_TBE_242.sh & ./script_scjack_540_TBE_243.sh & ./script_scjack_540_TBE_244.sh & ./script_scjack_540_TBE_245.sh & ./script_scjack_540_TBE_246.sh & ./script_scjack_540_TBE_247.sh & ./script_scjack_540_TBE_248.sh & ./script_scjack_540_TBE_249.sh & ./script_scjack_540_TBE_250.sh

echo "Commands 201-250 complete"

./script_scjack_540_TBE_251.sh & ./script_scjack_540_TBE_252.sh & ./script_scjack_540_TBE_253.sh & ./script_scjack_540_TBE_254.sh & ./script_scjack_540_TBE_255.sh & ./script_scjack_540_TBE_256.sh & ./script_scjack_540_TBE_257.sh & ./script_scjack_540_TBE_258.sh & ./script_scjack_540_TBE_259.sh & ./script_scjack_540_TBE_260.sh & ./script_scjack_540_TBE_261.sh & ./script_scjack_540_TBE_262.sh & ./script_scjack_540_TBE_263.sh & ./script_scjack_540_TBE_264.sh & ./script_scjack_540_TBE_265.sh & ./script_scjack_540_TBE_266.sh & ./script_scjack_540_TBE_267.sh & ./script_scjack_540_TBE_268.sh & ./script_scjack_540_TBE_269.sh & ./script_scjack_540_TBE_270.sh & ./script_scjack_540_TBE_271.sh & ./script_scjack_540_TBE_272.sh & ./script_scjack_540_TBE_273.sh & ./script_scjack_540_TBE_274.sh & ./script_scjack_540_TBE_275.sh & ./script_scjack_540_TBE_276.sh & ./script_scjack_540_TBE_277.sh & ./script_scjack_540_TBE_278.sh & ./script_scjack_540_TBE_279.sh & ./script_scjack_540_TBE_280.sh & ./script_scjack_540_TBE_281.sh & ./script_scjack_540_TBE_282.sh & ./script_scjack_540_TBE_283.sh & ./script_scjack_540_TBE_284.sh & ./script_scjack_540_TBE_285.sh & ./script_scjack_540_TBE_286.sh & ./script_scjack_540_TBE_287.sh & ./script_scjack_540_TBE_288.sh & ./script_scjack_540_TBE_289.sh & ./script_scjack_540_TBE_290.sh & ./script_scjack_540_TBE_291.sh & ./script_scjack_540_TBE_292.sh & ./script_scjack_540_TBE_293.sh & ./script_scjack_540_TBE_294.sh & ./script_scjack_540_TBE_295.sh & ./script_scjack_540_TBE_296.sh & ./script_scjack_540_TBE_297.sh & ./script_scjack_540_TBE_298.sh & ./script_scjack_540_TBE_299.sh & ./script_scjack_540_TBE_300.sh

echo "Commands 251-300 complete"

./script_scjack_540_TBE_301.sh & ./script_scjack_540_TBE_302.sh & ./script_scjack_540_TBE_303.sh & ./script_scjack_540_TBE_304.sh & ./script_scjack_540_TBE_305.sh & ./script_scjack_540_TBE_306.sh & ./script_scjack_540_TBE_307.sh & ./script_scjack_540_TBE_308.sh & ./script_scjack_540_TBE_309.sh & ./script_scjack_540_TBE_310.sh & ./script_scjack_540_TBE_311.sh & ./script_scjack_540_TBE_312.sh & ./script_scjack_540_TBE_313.sh & ./script_scjack_540_TBE_314.sh & ./script_scjack_540_TBE_315.sh & ./script_scjack_540_TBE_316.sh & ./script_scjack_540_TBE_317.sh & ./script_scjack_540_TBE_318.sh & ./script_scjack_540_TBE_319.sh & ./script_scjack_540_TBE_320.sh & ./script_scjack_540_TBE_321.sh & ./script_scjack_540_TBE_322.sh & ./script_scjack_540_TBE_323.sh & ./script_scjack_540_TBE_324.sh & ./script_scjack_540_TBE_325.sh & ./script_scjack_540_TBE_326.sh & ./script_scjack_540_TBE_327.sh & ./script_scjack_540_TBE_328.sh & ./script_scjack_540_TBE_329.sh & ./script_scjack_540_TBE_330.sh & ./script_scjack_540_TBE_331.sh & ./script_scjack_540_TBE_332.sh & ./script_scjack_540_TBE_333.sh & ./script_scjack_540_TBE_334.sh & ./script_scjack_540_TBE_335.sh & ./script_scjack_540_TBE_336.sh & ./script_scjack_540_TBE_337.sh & ./script_scjack_540_TBE_338.sh & ./script_scjack_540_TBE_339.sh & ./script_scjack_540_TBE_340.sh & ./script_scjack_540_TBE_341.sh & ./script_scjack_540_TBE_342.sh & ./script_scjack_540_TBE_343.sh & ./script_scjack_540_TBE_344.sh & ./script_scjack_540_TBE_345.sh & ./script_scjack_540_TBE_346.sh & ./script_scjack_540_TBE_347.sh & ./script_scjack_540_TBE_348.sh & ./script_scjack_540_TBE_349.sh & ./script_scjack_540_TBE_350.sh

echo "Commands 301-350 complete"

./script_scjack_540_TBE_351.sh & ./script_scjack_540_TBE_352.sh & ./script_scjack_540_TBE_353.sh & ./script_scjack_540_TBE_354.sh & ./script_scjack_540_TBE_355.sh & ./script_scjack_540_TBE_356.sh & ./script_scjack_540_TBE_357.sh & ./script_scjack_540_TBE_358.sh & ./script_scjack_540_TBE_359.sh & ./script_scjack_540_TBE_360.sh & ./script_scjack_540_TBE_361.sh & ./script_scjack_540_TBE_362.sh & ./script_scjack_540_TBE_363.sh & ./script_scjack_540_TBE_364.sh & ./script_scjack_540_TBE_365.sh & ./script_scjack_540_TBE_366.sh & ./script_scjack_540_TBE_367.sh & ./script_scjack_540_TBE_368.sh & ./script_scjack_540_TBE_369.sh & ./script_scjack_540_TBE_370.sh & ./script_scjack_540_TBE_371.sh & ./script_scjack_540_TBE_372.sh & ./script_scjack_540_TBE_373.sh & ./script_scjack_540_TBE_374.sh & ./script_scjack_540_TBE_375.sh & ./script_scjack_540_TBE_376.sh & ./script_scjack_540_TBE_377.sh & ./script_scjack_540_TBE_378.sh & ./script_scjack_540_TBE_379.sh & ./script_scjack_540_TBE_380.sh & ./script_scjack_540_TBE_381.sh & ./script_scjack_540_TBE_382.sh & ./script_scjack_540_TBE_383.sh & ./script_scjack_540_TBE_384.sh & ./script_scjack_540_TBE_385.sh & ./script_scjack_540_TBE_386.sh & ./script_scjack_540_TBE_387.sh & ./script_scjack_540_TBE_388.sh & ./script_scjack_540_TBE_389.sh & ./script_scjack_540_TBE_390.sh & ./script_scjack_540_TBE_391.sh & ./script_scjack_540_TBE_392.sh & ./script_scjack_540_TBE_393.sh & ./script_scjack_540_TBE_394.sh & ./script_scjack_540_TBE_395.sh & ./script_scjack_540_TBE_396.sh & ./script_scjack_540_TBE_397.sh & ./script_scjack_540_TBE_398.sh & ./script_scjack_540_TBE_399.sh & ./script_scjack_540_TBE_400.sh

echo "Commands 351-400 complete"

./script_scjack_540_TBE_401.sh & ./script_scjack_540_TBE_402.sh & ./script_scjack_540_TBE_403.sh & ./script_scjack_540_TBE_404.sh & ./script_scjack_540_TBE_405.sh & ./script_scjack_540_TBE_406.sh & ./script_scjack_540_TBE_407.sh & ./script_scjack_540_TBE_408.sh & ./script_scjack_540_TBE_409.sh & ./script_scjack_540_TBE_410.sh & ./script_scjack_540_TBE_411.sh & ./script_scjack_540_TBE_412.sh & ./script_scjack_540_TBE_413.sh & ./script_scjack_540_TBE_414.sh & ./script_scjack_540_TBE_415.sh & ./script_scjack_540_TBE_416.sh & ./script_scjack_540_TBE_417.sh & ./script_scjack_540_TBE_418.sh & ./script_scjack_540_TBE_419.sh & ./script_scjack_540_TBE_420.sh & ./script_scjack_540_TBE_421.sh & ./script_scjack_540_TBE_422.sh & ./script_scjack_540_TBE_423.sh & ./script_scjack_540_TBE_424.sh & ./script_scjack_540_TBE_425.sh & ./script_scjack_540_TBE_426.sh & ./script_scjack_540_TBE_427.sh & ./script_scjack_540_TBE_428.sh & ./script_scjack_540_TBE_429.sh & ./script_scjack_540_TBE_430.sh & ./script_scjack_540_TBE_431.sh & ./script_scjack_540_TBE_432.sh & ./script_scjack_540_TBE_433.sh & ./script_scjack_540_TBE_434.sh & ./script_scjack_540_TBE_435.sh & ./script_scjack_540_TBE_436.sh & ./script_scjack_540_TBE_437.sh & ./script_scjack_540_TBE_438.sh & ./script_scjack_540_TBE_439.sh & ./script_scjack_540_TBE_440.sh & ./script_scjack_540_TBE_441.sh & ./script_scjack_540_TBE_442.sh & ./script_scjack_540_TBE_443.sh & ./script_scjack_540_TBE_444.sh & ./script_scjack_540_TBE_445.sh & ./script_scjack_540_TBE_446.sh & ./script_scjack_540_TBE_447.sh & ./script_scjack_540_TBE_448.sh & ./script_scjack_540_TBE_449.sh & ./script_scjack_540_TBE_450.sh

echo "Commands 401-450 complete"

./script_scjack_540_TBE_451.sh & ./script_scjack_540_TBE_452.sh & ./script_scjack_540_TBE_453.sh & ./script_scjack_540_TBE_454.sh & ./script_scjack_540_TBE_455.sh & ./script_scjack_540_TBE_456.sh & ./script_scjack_540_TBE_457.sh & ./script_scjack_540_TBE_458.sh & ./script_scjack_540_TBE_459.sh & ./script_scjack_540_TBE_460.sh & ./script_scjack_540_TBE_461.sh & ./script_scjack_540_TBE_462.sh & ./script_scjack_540_TBE_463.sh & ./script_scjack_540_TBE_464.sh & ./script_scjack_540_TBE_465.sh & ./script_scjack_540_TBE_466.sh & ./script_scjack_540_TBE_467.sh & ./script_scjack_540_TBE_468.sh & ./script_scjack_540_TBE_469.sh & ./script_scjack_540_TBE_470.sh & ./script_scjack_540_TBE_471.sh & ./script_scjack_540_TBE_472.sh & ./script_scjack_540_TBE_473.sh & ./script_scjack_540_TBE_474.sh & ./script_scjack_540_TBE_475.sh & ./script_scjack_540_TBE_476.sh & ./script_scjack_540_TBE_477.sh & ./script_scjack_540_TBE_478.sh & ./script_scjack_540_TBE_479.sh & ./script_scjack_540_TBE_480.sh & ./script_scjack_540_TBE_481.sh & ./script_scjack_540_TBE_482.sh & ./script_scjack_540_TBE_483.sh & ./script_scjack_540_TBE_484.sh & ./script_scjack_540_TBE_485.sh & ./script_scjack_540_TBE_486.sh & ./script_scjack_540_TBE_487.sh & ./script_scjack_540_TBE_488.sh & ./script_scjack_540_TBE_489.sh & ./script_scjack_540_TBE_490.sh & ./script_scjack_540_TBE_491.sh & ./script_scjack_540_TBE_492.sh & ./script_scjack_540_TBE_493.sh & ./script_scjack_540_TBE_494.sh & ./script_scjack_540_TBE_495.sh & ./script_scjack_540_TBE_496.sh & ./script_scjack_540_TBE_497.sh & ./script_scjack_540_TBE_498.sh & ./script_scjack_540_TBE_499.sh & ./script_scjack_540_TBE_500.sh

echo "Commands 451-500 complete"


echo "ALL COMMANDS COMPLETE!!!!!"

exit