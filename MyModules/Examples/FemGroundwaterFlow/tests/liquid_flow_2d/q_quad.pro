<?xml version="1.0"?>
<ogs6>
<coupling>
	<P algorithm="Serial" convergence="FemFunctionConvergenceCheck" max_itr="1" epsilon="1e-4">
		<out>PRESSURE1</out>
		<out>VELOCITY</out>
		<problems>
			<M type="LIQUID_FLOW">
				<out>PRESSURE1</out>
			</M>
			<M type="PRESSURE_TO_ELEMENT_VELOCITY">
				<in>PRESSURE1</in>
				<out>VELOCITY</out>
			</M>
		</problems>
	</P>
</coupling>
</ogs6>

