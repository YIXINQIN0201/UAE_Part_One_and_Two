<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis styleCategories="Symbology" version="3.44.7-Solothurn">
  <pipe-data-defined-properties>
    <Option type="Map">
      <Option type="QString" value="" name="name"/>
      <Option name="properties"/>
      <Option type="QString" value="collection" name="type"/>
    </Option>
  </pipe-data-defined-properties>
  <pipe>
    <provider>
      <resampling enabled="false" zoomedInResamplingMethod="nearestNeighbour" maxOversampling="2" zoomedOutResamplingMethod="nearestNeighbour"/>
    </provider>
    <rasterrenderer alphaBand="-1" nodataColor="" type="paletted" opacity="0.739" band="1">
      <rasterTransparency>
        <singleValuePixelList>
          <pixelListEntry percentTransparent="0" min="nan" max="nan"/>
        </singleValuePixelList>
      </rasterTransparency>
      <minMaxOrigin>
        <limits>None</limits>
        <extent>WholeRaster</extent>
        <statAccuracy>Estimated</statAccuracy>
        <cumulativeCutLower>0.02</cumulativeCutLower>
        <cumulativeCutUpper>0.98</cumulativeCutUpper>
        <stdDevFactor>2</stdDevFactor>
      </minMaxOrigin>
      <colorPalette>
        <paletteEntry label="Degarding" alpha="255" value="1" color="#c46c61"/>
      </colorPalette>
      <colorramp type="randomcolors" name="[source]">
        <Option/>
      </colorramp>
    </rasterrenderer>
    <brightnesscontrast brightness="0" gamma="1" contrast="0"/>
    <huesaturation colorizeBlue="97" saturation="0" grayscaleMode="0" colorizeRed="196" invertColors="0" colorizeGreen="108" colorizeStrength="100" colorizeOn="0"/>
    <rasterresampler maxOversampling="2"/>
    <resamplingStage>resamplingFilter</resamplingStage>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
