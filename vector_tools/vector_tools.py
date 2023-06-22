# coding: utf-8

"""
DOCME

***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

import fiona
import fiona.crs
from shapely.geometry import shape, box
from shapely.ops import unary_union
from rtree import index
from shapely.geometry import LineString, MultiLineString, mapping, Point
import numpy as np

def FeatureInPolygonWithDistance(linestring_path, polygon_path, output_path, distance_threshold):

    selected_features = []

    # Create spatial index for polygon layer
    polygon_index = index.Index()
    with fiona.open(polygon_path, 'r') as polygon_layer:
        for idx, polygon_feature in polygon_layer.items():
            polygon = shape(polygon_feature['geometry'])
            polygon_index.insert(idx, polygon.bounds)

    # Create spatial index for linestring layer
    linestring_index = index.Index()
    with fiona.open(linestring_path, 'r') as linestring_layer:
        for idx, linestring_feature in linestring_layer.items():
            linestring = shape(linestring_feature['geometry'])
            linestring_index.insert(idx, linestring.bounds)

    with fiona.open(polygon_path, 'r') as polygon_layer, fiona.open(linestring_path, 'r') as linestring_layer:
        options = dict(
                driver=linestring_layer.driver,
                schema=linestring_layer.schema.copy(),
                crs=linestring_layer.crs)

        for polygon_feature in polygon_layer:
            # print(polygon_feature)
            polygon = shape(polygon_feature['geometry'])

            # Get potential linestrings within the bounding box of the polygon
            potential_linestrings = list(linestring_index.intersection(polygon.bounds))

            for linestring_id  in potential_linestrings:
                linestring_feature = linestring_layer[linestring_id]
                # print(linestring_feature)
                linestring = shape(linestring_feature['geometry'])

                if linestring.intersects(polygon):
                    clipped_linestring = linestring.intersection(polygon)
                    linestring_length = clipped_linestring.length
                    if linestring_length > distance_threshold:
                        selected_features.append(linestring_feature)
    
    # Create a new GeoPackage file and write the selected features to it
    with fiona.open(output_path, 'w', **options) as output_layer:
        for feature in selected_features:
            output_layer.write(feature)


def convert_multi_line_to_line(input_file, output_file):
    # Open the input file with Fiona
    with fiona.open(input_file, 'r') as src:
        # Create a new output file with the desired schema
        schema=src.schema.copy()
        options = dict(
                driver=src.driver,
                schema=schema,
                crs=src.crs)
        with fiona.open(output_file, 'w', **options) as dst:
            # Iterate over the features in the input file
            for feature in src:
                # Get the geometry of the feature
                geometry = shape(feature['geometry'])
                
                # If it's a MultiLineString, convert it to LineString(s)
                if geometry.geom_type == 'MultiLineString':
                    for line in geometry:
                        # Create a new feature with the LineString geometry
                        new_feature = feature.copy()
                        new_feature['geometry'] = mapping(line)
                        
                        # Write the new feature to the output file
                        dst.write(new_feature)
                
                # If it's already a LineString, write it directly to the output file
                elif geometry.geom_type == 'LineString':
                    dst.write(feature)


def identify_network_nodes(input_file, output_file, nodes_file):
    with fiona.open(input_file, 'r') as source:
        driver = source.driver
        source_crs = source.crs
        schema = source.schema.copy()

        nodes_schema = {
            'geometry': 'Point',
            'properties': {'GID': 'int:10'}
        }

        output_schema = schema.copy()
        output_schema['properties']['NODEA'] = 'int:10'
        output_schema['properties']['NODEB'] = 'int:10'

        with fiona.open(output_file, 'w', driver, output_schema, crs=source_crs) as output, \
             fiona.open(nodes_file, 'w', driver, nodes_schema, crs=source_crs) as nodes:
            coordinates = []

            for feature in source:
                line = LineString(feature['geometry']['coordinates'])
                a = line.coords[0]
                b = line.coords[-1]
                coordinates.append(Point(a))
                coordinates.append(Point(b))

            # Quantization
            minx = min(point.x for point in coordinates)
            miny = min(point.y for point in coordinates)
            maxx = max(point.x for point in coordinates)
            maxy = max(point.y for point in coordinates)

            kx = 1 if minx == maxx else maxx - minx
            ky = 1 if miny == maxy else maxy - miny
            sx = kx / 1e8
            sy = ky / 1e8

            coordinates = [Point((point.x - minx) / sx, (point.y - miny) / sy) for point in coordinates]

            # Build Endpoints Index
            point_index = {}
            gid = 0

            for coordinate in coordinates:
                if coordinate not in point_index:
                    point_index[coordinate] = gid
                    point = Point(coordinate.x * sx + minx, coordinate.y * sy + miny)
                    nodes.write({
                        'geometry': {'type': 'Point', 'coordinates': (point.x, point.y)},
                        'properties': {'GID': gid}
                    })
                    gid += 1

            # Output Lines with Node Attributes
            for feature in source:
                geometry = LineString(feature['geometry']['coordinates'])

                if geometry.is_simple:
                    a = geometry.coords[0]
                    b = geometry.coords[-1]
                    node_a = point_index[Point((a[0] - minx) / sx, (a[1] - miny) / sy)]
                    node_b = point_index[Point((b[0] - minx) / sx, (b[1] - miny) / sy)]
                    feature['properties']['NODEA'] = node_a
                    feature['properties']['NODEB'] = node_b
                    output.write(feature)
                else:
                    for part in geometry:
                        a = part.coords[0]
                        b = part.coords[-1]
                        node_a = point_index[Point((a[0] - minx) / sx, (a[1] - miny) / sy)]
                        node_b = point_index[Point((b[0] - minx) / sx, (b[1] - miny) / sy)]
                        output.write({
                            'geometry': {'type': 'LineString', 'coordinates': list(part.coords)},
                            'properties': {'NODEA': node_a, 'NODEB': node_b}
                        })

def ExtractBylocation(input_file, mask_file, output_file, method):
    selected_features = []

    # Create spatial index for input layer
    input_index = index.Index()
    with fiona.open(input_file, 'r') as input_layer:
        for idx, input_feature in input_layer.items():
            in_feat = shape(input_feature['geometry'])
            input_index.insert(idx, in_feat.bounds)

    # Create spatial index for mask layer
    mask_index = index.Index()
    with fiona.open(mask_file, 'r') as mask_layer:
        for idx, mask_feature in mask_layer.items():
            mask_feat = shape(mask_feature['geometry'])
            mask_index.insert(idx, mask_feat.bounds)
    
    with fiona.open(mask_file, 'r') as mask_layer, fiona.open(input_file, 'r') as input_layer:
        options = dict(
                driver=input_layer.driver,
                schema=input_layer.schema.copy(),
                crs=input_layer.crs)

        for mask_feature in mask_layer:
            # print(polygon_feature)
            mask = shape(mask_feature['geometry'])

            # Get potential linestrings within the bounding box of the polygon
            potential_input = list(input_index.intersection(mask.bounds))

            for input_id  in potential_input:
                input_feature = input_layer[input_id]
                # print(linestring_feature)
                input = shape(input_feature['geometry'])

                if method == 'intersects':
                    if mask.intersects(input):
                        selected_features.append(input_feature)
                if method == 'contains':
                    if mask.contains(input):
                        selected_features.append(input_feature)
    # Create a new GeoPackage file and write the selected features to it
    with fiona.open(output_file, 'w', **options) as output_layer:
        for feature in selected_features:
            output_layer.write(feature)

def ExtractByBoundMask(input_file, mask_file, output_file):
    selected_features = []

    # # Create spatial index for input layer
    # input_index = index.Index()
    # with fiona.open(input_file, 'r') as input_layer:
    #     for idx, input_feature in input_layer.items():
    #         in_feat = shape(input_feature['geometry'])
    #         input_index.insert(idx, in_feat.bounds)
    
    with fiona.open(mask_file, 'r') as mask_layer:
        
        polygons = [shape(feature["geometry"]) for feature in mask_layer]

        # Merge the polygons into a single polygon
        merged_polygon = unary_union(polygons)
        # Calculate the bounding box of the merged polygon
        bounds = box(*merged_polygon.bounds)

        with fiona.open(input_file, 'r') as input_layer:
            options = dict(
                driver=input_layer.driver,
                schema=input_layer.schema.copy(),
                crs=input_layer.crs)
            for feat in input_layer:
                geom = shape(feat["geometry"])
                if geom.intersects(bounds):
                    selected_features.append(feat)
        # Create a new GeoPackage file and write the selected features to it
        with fiona.open(output_file, 'w', **options) as output_layer:
            for select in selected_features:
                output_layer.write(select) 

# work from python-fct not finish (for memory)
def JoinNetworkAttributes(params, tileset = 'default'):

    sources_confluences = params.sources_confluences.filename(tileset=None)
    network_identified_strahler = params.network_identified_strahler.filename(tileset=tileset)
    rhts = params.rhts.filename(tileset=tileset)

    strahler_field_name = 'strahler'
    distance_threshold = 50

    # get the fields in the schema from sources_confluences not existing in network_identified_strahler 
    with fiona.open(sources_confluences, 'r') as source:
        with fiona.open(network_identified_strahler, 'r') as network:
            schema_source_prop = source.schema['properties']
            schema_network_prop = network.schema['properties']
            new_schema_network = network.schema.copy()
            for key, value in schema_source_prop.items():
                if key not in schema_network_prop:
                    new_schema_network['properties'][key] = value

    # get sources_confluences with strahler == 1
    with fiona.open(sources_confluences, 'r') as source:
        strahler1 = [feature for feature in source if feature['properties'][strahler_field_name] == 1]


    # Create an R-tree index
    strahler1_index = index.Index()
    # populate index
    for i, point in enumerate(strahler1):
        geom = shape(point['geometry'])
        strahler1_index.insert(i, geom.bounds)


    with fiona.open(network_identified_strahler, 'r') as network:
        
        options = dict(
            schema = new_schema_network,
            driver = network.driver,
            crs = network.crs
        )
        with fiona.open(rhts, 'w', **options) as output:
        
            for line in network:
                # initialisation to get nearest_point
                nearest_point = None
                nearest_distance = float('inf')
                line_properties = line['properties']
                # get lione with strahler == 1
                if line['properties'][strahler_field_name] == 1:
                    geometry = shape(line['geometry'])
                    # get line first point coordinates
                    first_point = Point(geometry.coords[0])
                    # create buffer around first point
                    point_buffer = first_point.buffer(distance_threshold)
                    # get all the point from source that intersect with buffer with index
                    potential_matches = [idx for idx in strahler1_index.intersection(point_buffer.bounds)] # store idx of the index
                    
                    if potential_matches:
                        # seach with index in potential_matches
                        for idx in potential_matches:
                            source_potential_geom = shape(strahler1[idx]['geometry'])
                            # iterate through potential_matches strahler1, calculate distance and get the shortest by updating nearest_distance if the current strahler1 is closer
                            if first_point.distance(source_potential_geom) < nearest_distance:
                                nearest_point = strahler1[idx]
                                nearest_distance = first_point.distance(source_potential_geom)
                        nearest_point_properties = nearest_point['properties']

                        # get the fields in the properties from sources_confluences not existing in network_identified_strahler 
                        output_properties = line_properties
                        for key, value in nearest_point_properties.items():
                            print(key)
                            if key not in line_properties:
                                output_properties[key] = value
                        print(output_properties)
                        rhts_feature = {
                                        'type': 'Feature',
                                        'properties': output_properties,
                                        'geometry': line['geometry'],
                                    }
                        output.write(rhts_feature)