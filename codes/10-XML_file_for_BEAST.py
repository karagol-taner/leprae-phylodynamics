import os

fasta_in = "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/alignment_fixed.fasta"
xml_out = "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/lep_beast_full_v27.xml"

print("--- Generating FULL-SCALE BEAST v2.7.7 XML (RandomTree FIX) ---")

with open(fasta_in, 'r') as f:
    lines = f.readlines()

taxa = []
for i in range(0, len(lines), 2):
    name = lines[i].strip()[1:]
    seq = lines[i+1].strip()
    taxa.append((name, seq))

log_name = xml_out.replace(".xml", ".log")
trees_name = xml_out.replace(".xml", ".trees")

trait_entries = []
for name, _ in taxa:
    if "SRR847" in name:
        year = "1050"
    elif "SRR367" in name:
        year = "2016"
    else:
        year = "2020"
    trait_entries.append(f"{name}={year}")
trait_string = ",\n                    ".join(trait_entries)

with open(xml_out, 'w') as x:
    x.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
    x.write('<beast beautitemplate="Standard" namespace="beast.pkgmgmt:beast.base.core:beast.base.inference:beast.base.inference.operator:beast.base.inference.distribution:beast.base.evolution.alignment:beast.base.evolution.tree.coalescent:beast.base.util:beast.base.evolution.nuc:beast.base.evolution.operator:beast.base.evolution.sitemodel:beast.base.evolution.substitutionmodel:beast.base.evolution.likelihood:beast.base.evolution.tree:beast.base.evolution.speciation:beast.base.evolution.branchratemodel" version="2.7">\n')

    # 1. Alignment
    x.write('    <data id="alignment" name="alignment" dataType="nucleotide">\n')
    for name, seq in taxa:
        x.write(f'        <sequence id="seq_{name}" taxon="{name}" totalcount="4" value="{seq}"/>\n')
    x.write('    </data>\n')

    # 2. MCMC
    x.write('    <run id="mcmc" spec="beast.base.inference.MCMC" chainLength="50000000">\n')
    x.write('        <state id="state" storeEvery="10000">\n')

    x.write('            <tree id="Tree.t:alignment" name="stateNode" spec="beast.base.evolution.tree.Tree">\n')
    x.write(f'                <trait id="dateTrait.t:alignment" spec="beast.base.evolution.tree.TraitSet" traitname="date" value="{trait_string}">\n')
    x.write('                    <taxa id="TaxonSet.alignment" spec="beast.base.evolution.alignment.TaxonSet">\n')
    x.write('                        <alignment idref="alignment"/>\n')
    x.write('                    </taxa>\n')
    x.write('                </trait>\n')
    x.write('            </tree>\n')

    x.write('            <parameter id="birthRate.t:alignment" name="stateNode" lower="0.0">1.0</parameter>\n')
    x.write('            <parameter id="clockRate.c:alignment" name="stateNode" lower="1.0E-10" upper="1.0E-5">1.0E-8</parameter>\n')

    for rate in ['AC', 'AG', 'AT', 'CG', 'CT', 'GT']:
        x.write(f'            <parameter id="rate{rate}.s:alignment" name="stateNode" lower="0.0">1.0</parameter>\n')
    x.write('            <parameter id="freqParameter.s:alignment" dimension="4" name="stateNode" lower="0.0" upper="1.0">0.25 0.25 0.25 0.25</parameter>\n')
    x.write('        </state>\n')

    # === RandomTree FIX ===
    # Changed spec to beast.base.evolution.tree.coalescent.RandomTree
    x.write('        <init id="RandomTree.t:alignment" spec="beast.base.evolution.tree.coalescent.RandomTree" estimate="false" initial="@Tree.t:alignment" taxa="@alignment">\n')
    x.write('            <populationModel id="ConstantPopulation0.t:alignment" spec="beast.base.evolution.tree.coalescent.ConstantPopulation">\n')
    x.write('                <parameter id="randomPopSize.t:alignment" name="popSize">1.0</parameter>\n')
    x.write('            </populationModel>\n')
    x.write('        </init>\n')

    # 3. Posterior
    x.write('        <distribution id="posterior" spec="beast.base.inference.CompoundDistribution">\n')
    x.write('            <distribution id="likelihood" spec="beast.base.inference.CompoundDistribution" useThreads="true">\n')
    x.write('                <distribution id="treeLikelihood.alignment" spec="beast.base.evolution.likelihood.ThreadedTreeLikelihood">\n')
    x.write('                    <data id="filtered_alignment" spec="beast.base.evolution.alignment.FilteredAlignment" data="@alignment" constantSiteWeights="754321 882145 881234 752345" filter="1-"/>\n')
    x.write('                    <tree idref="Tree.t:alignment"/>\n')
    x.write('                    <siteModel id="SiteModel.s:alignment" spec="beast.base.evolution.sitemodel.SiteModel">\n')
    x.write('                        <substModel id="gtr.s:alignment" spec="beast.base.evolution.substitutionmodel.GTR">\n')
    for rate in ['AC', 'AG', 'AT', 'CG', 'CT', 'GT']:
        x.write(f'                            <parameter idref="rate{rate}.s:alignment" name="rate{rate}"/>\n')
    x.write('                            <frequencies id="estimatedFreqs.s:alignment" spec="beast.base.evolution.substitutionmodel.Frequencies">\n')
    x.write('                                <parameter idref="freqParameter.s:alignment" name="frequencies"/>\n')
    x.write('                            </frequencies>\n')
    x.write('                        </substModel>\n')
    x.write('                    </siteModel>\n')
    x.write('                    <branchRateModel id="StrictClock.c:alignment" spec="beast.base.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.c:alignment"/>\n')
    x.write('                </distribution>\n')
    x.write('            </distribution>\n')

    # 3b. Priors
    x.write('            <distribution id="prior" spec="beast.base.inference.CompoundDistribution">\n')
    x.write('                <distribution id="YuleModel.t:alignment" spec="beast.base.evolution.speciation.YuleModel" birthDiffRate="@birthRate.t:alignment" tree="@Tree.t:alignment"/>\n')
    x.write('                <distribution id="ClockPrior.c:alignment" spec="beast.base.inference.distribution.Prior" name="distribution" x="@clockRate.c:alignment">\n')
    x.write('                    <distr id="Uniform.clock" spec="beast.base.inference.distribution.Uniform" lower="1.0E-10" upper="1.0E-5"/>\n')
    x.write('                </distribution>\n')
    x.write('                <distribution id="YuleBirthRatePrior.t:alignment" spec="beast.base.inference.distribution.Prior" name="distribution" x="@birthRate.t:alignment">\n')
    x.write('                    <distr id="LogNormal.birthRate" spec="beast.base.inference.distribution.LogNormalDistributionModel" M="1.0" S="1.25"/>\n')
    x.write('                </distribution>\n')

    for rate in ['AC', 'AG', 'AT', 'CG', 'CT', 'GT']:
        x.write(f'                <distribution id="rate{rate}Prior.s:alignment" spec="beast.base.inference.distribution.Prior" name="distribution" x="@rate{rate}.s:alignment">\n')
        x.write(f'                    <distr id="Gamma.rate{rate}" spec="beast.base.inference.distribution.Gamma" alpha="0.05" beta="10.0"/>\n')
        x.write(f'                </distribution>\n')

    x.write('                <distribution id="FrequenciesPrior.s:alignment" spec="beast.base.inference.distribution.Prior" name="distribution" x="@freqParameter.s:alignment">\n')
    x.write('                    <distr id="Dirichlet.freqs" spec="beast.base.inference.distribution.Dirichlet" alpha="4.0 4.0 4.0 4.0"/>\n')
    x.write('                </distribution>\n')
    x.write('            </distribution>\n')
    x.write('        </distribution>\n')

    # 4. Operators
    x.write('        <operator id="treeScaler.t:alignment" spec="beast.base.evolution.operator.ScaleOperator" scaleFactor="0.5" tree="@Tree.t:alignment" weight="3.0"/>\n')
    x.write('        <operator id="treeRootScaler.t:alignment" spec="beast.base.evolution.operator.ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:alignment" weight="3.0"/>\n')
    x.write('        <operator id="UniformOperator.t:alignment" spec="beast.base.evolution.operator.Uniform" tree="@Tree.t:alignment" weight="30.0"/>\n')
    x.write('        <operator id="SubtreeSlide.t:alignment" spec="beast.base.evolution.operator.SubtreeSlide" tree="@Tree.t:alignment" weight="15.0"/>\n')
    x.write('        <operator id="narrow.t:alignment" spec="beast.base.evolution.operator.Exchange" tree="@Tree.t:alignment" weight="15.0"/>\n')
    x.write('        <operator id="wide.t:alignment" spec="beast.base.evolution.operator.Exchange" isNarrow="false" tree="@Tree.t:alignment" weight="3.0"/>\n')
    x.write('        <operator id="WilsonBalding.t:alignment" spec="beast.base.evolution.operator.WilsonBalding" tree="@Tree.t:alignment" weight="3.0"/>\n')
    x.write('        <operator id="birthRateScaler.t:alignment" spec="beast.base.evolution.operator.ScaleOperator" parameter="@birthRate.t:alignment" scaleFactor="0.75" weight="3.0"/>\n')
    x.write('        <operator id="StrictClockRateScaler.c:alignment" spec="beast.base.evolution.operator.ScaleOperator" parameter="@clockRate.c:alignment" scaleFactor="0.75" weight="3.0"/>\n')

    x.write('        <operator id="strictClockUpDownOperator.c:alignment" spec="beast.base.inference.operator.UpDownOperator" scaleFactor="0.75" weight="3.0">\n')
    x.write('            <up idref="clockRate.c:alignment"/>\n')
    x.write('            <down idref="Tree.t:alignment"/>\n')
    x.write('        </operator>\n')

    for i, rate in enumerate(['AC', 'AG', 'AT', 'CG', 'CT', 'GT'], 1):
        x.write(f'        <operator id="GTRRatesScaler{i}.s:alignment" spec="beast.base.evolution.operator.ScaleOperator" parameter="@rate{rate}.s:alignment" scaleFactor="0.5" weight="0.1"/>\n')
    x.write('        <operator id="FrequenciesExchanger.s:alignment" spec="beast.base.inference.operator.DeltaExchangeOperator" parameter="@freqParameter.s:alignment" delta="0.01" weight="0.1"/>\n')

    # 5. Loggers
    x.write(f'        <logger id="tracelog" spec="beast.base.inference.Logger" fileName="{log_name}" logEvery="1000">\n')
    x.write('            <log idref="posterior"/>\n')
    x.write('            <log idref="likelihood"/>\n')
    x.write('            <log idref="prior"/>\n')
    x.write('            <log id="TreeHeight.t:alignment" spec="beast.base.evolution.tree.TreeHeightLogger" tree="@Tree.t:alignment"/>\n')
    x.write('            <log idref="birthRate.t:alignment"/>\n')
    x.write('            <log idref="clockRate.c:alignment"/>\n')
    for rate in ['AC', 'AG', 'AT', 'CG', 'CT', 'GT']:
        x.write(f'            <log idref="rate{rate}.s:alignment"/>\n')
    x.write('            <log idref="freqParameter.s:alignment"/>\n')
    x.write('        </logger>\n')
    x.write('        <logger id="screenlog" spec="beast.base.inference.Logger" logEvery="1000">\n')
    x.write('            <log idref="posterior"/>\n')
    x.write('            <log idref="likelihood"/>\n')
    x.write('            <log idref="prior"/>\n')
    x.write('        </logger>\n')
    x.write(f'        <logger id="treelog" spec="beast.base.inference.Logger" fileName="{trees_name}" logEvery="10000" mode="tree">\n')
    x.write('            <log idref="Tree.t:alignment"/>\n')
    x.write('        </logger>\n')

    x.write('    </run>\n')
    x.write('</beast>\n')

print(f"✅ XML generated successfully: {xml_out}")
