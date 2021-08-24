### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ c25b5f6c-61a4-11eb-2c34-cbf860b5d314
using DataFrames, Dates, Dierckx, Distributions, GeoJSON, GeoInterface, Glob, HTTP, JSON, NCDatasets, PlutoUI, PolygonOps, StatsPlots, Statistics

# ╔═╡ 6659c6d4-7182-11eb-3b65-d939a40b4325
html"<button onclick=present()>Present</button>"

# ╔═╡ 1c2705ee-61a4-11eb-33f9-835cc5caa8a8
md"# TROPOMI/OCO-2 SIF DEMO"

# ╔═╡ b7921386-7186-11eb-20d2-d973f7940a35
md" **In this notebook we will:**
> - **Read L2 satellite data of solar-induced chlorophyll fluorescence (SIF) inferred from TROPOMI and OCO-2**

> - **Create time series over spatial features**

> - **Generate spatial composites**

> - **Evaluate the measurement uncertainties**"

# ╔═╡ e7799730-7187-11eb-02a8-e7f76ed6081c
md"## Start-up Julia"

# ╔═╡ b434a196-61a4-11eb-19fe-976e23dec488
md"**Loading required packages** which provide necessary functionality (e.g., reading certain data types)"

# ╔═╡ a4854224-7184-11eb-291a-e3667b3c33e2
md"**Switch to Plotly backend** to make plots interactive:"

# ╔═╡ 5250d5d2-718f-11eb-03ab-5d29ef4ce6ad
plotly()

# ╔═╡ 34be05d8-7185-11eb-3bd3-d1e24885ed75
md"## Download data"

# ╔═╡ 8cc541f6-7185-11eb-3034-675a8bf6ff7c
md"
> **TROPOMI data can be downloaded from our [ftp server](ftp://fluo.gps.caltech.edu/data/tropomi)** "

# ╔═╡ 2be873b4-718d-11eb-1a33-a742dd4610ff
LocalResource("./src/images/data_calendar.png", :width => 400)

# ╔═╡ bb5cb8a2-718d-11eb-3694-cbe3adf78a8d
md"
> **OCO-2 data can be downloaded from our [ftp server](ftp://fluo.gps.caltech.edu/data/OCO2/sif_lite_B8100/) or [GES DISC](https://disc.gsfc.nasa.gov/datasets/OCO2_L2_Lite_SIF_10r/summary?keywords=oco2%20sif%20lite)** "

# ╔═╡ a33e9ac2-718c-11eb-18d4-b90c627f52ec
md"## After Downloading the data"

# ╔═╡ 8f48bb84-718e-11eb-333e-8f1f2df863d9
md" Satellite data is often shared via .NC (netCDF, network Common Data Form) files. OCO-2/3 and TROPOMI data is stored in daily files, while TROPOMI files on our ftp server are bundled to monthly *.tar.gz archives. Unpacking the data is in this case the first step: ```tar -xvf *.tar.gz```
"

# ╔═╡ 38f5b8de-7190-11eb-0437-497aa85610f0
md"> **Satellite data sets have typically a large volume. In case of TROPOMI, one day includes about 4M observations over land with a size of up to 150MB (one month ~3.5GB).**"

# ╔═╡ 2df3dee4-61a4-11eb-2003-5b56cfb49d01
md"**Define location of the downloaded and unpacked files:**"

# ╔═╡ 5b5b38b4-61a4-11eb-284b-e329267e0fb2
## TROPOMI:
#sifDataDir = "/net/fluo/data2/data/TROPOMI_SIF740nm/original"
## OCO-2:
sifDataDir = "C:/Users/Owner/Downloads/NASA"

# ╔═╡ 4f00b9e0-61a4-11eb-0d33-f3ec82cd5396
md"**Get a list of available files:**"

# ╔═╡ 2f0d9a3a-61a5-11eb-2fc1-650e465800c8
## TROPOMI:
#allFiles = glob("*.nc",sifDataDir)
## OCO-2:
allFiles = glob("*.nc4",sifDataDir)

# ╔═╡ 62fa453c-61a5-11eb-0b39-015a9ff4155c
md" **Convert filenames to dates:**"

# ╔═╡ a18e5c16-61c3-11eb-37cf-edd1266f9dde
begin
	getTROPOMIdates = function(allFiles)
		## TROPOMI:
		FilePattern = "TROPO_SIF_YYYY-MM-DD_*.nc"	
		posY = findfirst("YYYY",FilePattern)
		posM = findfirst("MM",FilePattern)
		posD = findfirst("DD",FilePattern)
		allDates=Vector{Date}()
		for file in basename.(allFiles)
			push!(allDates,Date(file[posY]*file[posM]*file[posD], "yyyymmdd"))
		end
		return(allDates)
	end

	getOCO2dates = function(allFiles)
		FilePattern = "oco2_LtSIF_YYMMDD*.nc4"
		posY = findfirst("YY",FilePattern)
		posM = findfirst("MM",FilePattern)
		posD = findfirst("DD",FilePattern)
	
		allDates=Vector{Date}()
		for file in basename.(allFiles)
			push!(allDates,Date("20"*file[posY]*file[posM]*file[posD], "yyyymmdd"))
		end
		return(allDates)
	end


end

# ╔═╡ 179972c6-7257-11eb-3a3d-45c672612fb3
allDates = getOCO2dates(allFiles)

# ╔═╡ bc20f27e-61d1-11eb-2497-179d53ab2360
md" **Select date range:**"

# ╔═╡ dfaee354-61d1-11eb-2a88-8bb2002c2136
begin
	#startDate = Date(2020,1,1)
	#stopDate  = Date(2020,1,15)
	startDate = Date(2020,12,30)
	stopDate  = Date(2021,01,30)
	FileInd = findall(x -> x >= startDate && x <= stopDate, allDates)
end

# ╔═╡ a8dc819e-7188-11eb-2c1b-8998b48d99f5
md"## Create and read spatial features"

# ╔═╡ cfd6f918-7189-11eb-33f4-c5133381580b
md" **Defining Regions of interest:**

> One simple way to search for interesting features is using [Google Earth](https://www.google.com/earth/).
"

# ╔═╡ c475fc80-718a-11eb-2c79-3f45c09aaca4
LocalResource("./src/images/google_screenshot.png")

# ╔═╡ 2680b0a4-718c-11eb-2bad-137068539830
md"## Converting Google's KML files"

# ╔═╡ 471fe14a-718c-11eb-23a2-a3c5ed97e288
md" The Keyhole Markup Language (**KML**) uses a tag-based structure with nested elements and attributes consistent with the Extensible Markup Language (**XML**)."

# ╔═╡ a748c364-71a2-11eb-21fb-2ddc00692957
LocalResource("./src/images/kml_screenshot.png")

# ╔═╡ a8564f6a-71a2-11eb-1797-67abf06ac8c1
md"
Markup languages generally define a set of rules for encoding documents in a format that is both human-readable and machine-readable. Another widely used format with this property is the JavaScript Object Notation (**JSON**) and **GeoJSON**, an open standard for encoding geographic data structures."


# ╔═╡ d753c838-71a2-11eb-151d-6394549ab869
LocalResource("./src/images/json_screenshot.png")

# ╔═╡ d7f995ce-71a2-11eb-3e35-b113a803a9fb
md"
There are plenty of tools like this [one](https://mygeodata.cloud/converter/kml-to-geojson) to convert **KML** to **GeoJSON** files, which we will use in the following."

# ╔═╡ 91b73620-7193-11eb-070f-573411f9e611
md"## Reading GeoJSON files"

# ╔═╡ df1138ee-61d5-11eb-34b6-67ce049efb74
testROI=GeoJSON.read(read("C:/Users/Owner/Downloads/OCO-2.geojson"))

# ╔═╡ 14e5f60c-6687-11eb-37f6-fbf0d9061530
testROI.features[1].properties["name"]

# ╔═╡ ad814292-65d9-11eb-1390-35f78905835c
### Function to get coordinates w/o altitude for plotting and to make use of "inpolygon" from PolygonOps.jl later on:
getLonLat = function(x::Feature)
	LonLat = []
	for i in 1:length(x.geometry.coordinates[1])
		LonLat = push!(LonLat,x.geometry.coordinates[1][i][1:2])
	end
	return(Tuple.(LonLat))
end

# ╔═╡ e083f3fc-65da-11eb-290c-8df8464a4ff1
begin
	vegPolygon = getLonLat(testROI.features[1])
	Plots.plot(vegPolygon, legend=false)
	vegPolygon = getLonLat(testROI.features[1])
	Plots.plot(vegPolygon, legend=false)
	#Plots.plot!(mbr)
end

# ╔═╡ 53d0a158-62a6-11eb-3f8b-af81223169ec
md"## Compute the minimum boundary box (mbr)"

# ╔═╡ 213784d6-7198-11eb-16a2-254812fbed1d
md" > This step is necessary to narrow down the huge data volume."

# ╔═╡ b88cdc64-6509-11eb-24c7-fd78be2de33c
# accessing the coordinates works like this:
size(testROI.features[1].geometry.coordinates[1])

# ╔═╡ 80cc3a5c-629f-11eb-3807-f51439287448
getMBReps = function(x::FeatureCollection, epsLatLon=0.2)
	# MBR: minimum boundary box
	box   = zeros(4)
 	for ifeat in 1:length(x.features)
	
		tmp = x.features[ifeat].geometry.coordinates[1] |> DataFrame([126.17,33.2741],[126.263,33.1973],[126.285,33.1912],[126.323,33.2365],[126.601,33.2347],[126.841,33.3098],[126.943,33.4557],[126.902,33.5254],[126.778,33.5683],[126.316,33.4708],[126.165,33.3386],[126.17,33.2741], :auto)
		
		if ifeat==1
			box[1] = minimum(tmp[1,:])
			box[2] = maximum(tmp[1,:])
			box[3] = minimum(tmp[2,:])
			box[4] = maximum(tmp[2,:])
		end
		
		if ifeat > 1
			if minimum(tmp[1,:]) < box[1]
				box[1] = minimum(tmp[1,:])
			end
			if minimum(tmp[2,:]) < box[3] 
				box[3] = minimum(tmp[2,:])
			end
			if maximum(tmp[1,:]) > box[2] 
				box[2] = maximum(tmp[1,:])
			end
			if maximum(tmp[2,:]) > box[4] 
				box[4] = maximum(tmp[2,:])
			end
		end
	
	end
	box[1] = box[1]-epsLatLon
	box[2] = box[2]+epsLatLon
	box[3] = box[3]-epsLatLon
	box[4] = box[4]+epsLatLon

	return(box)
end

# ╔═╡ b1b2bb48-64fd-11eb-1ce8-274175a2d8ff
mbr=getMBReps(testROI)

# ╔═╡ b3fa469c-61d4-11eb-17ca-31a9319d9560
md"## Defining which parts of the NC files are relevant"

# ╔═╡ 50938100-62af-11eb-0027-2bc902261db2
# look at what's available in the NC files:
let
	ds = Dataset(allFiles[1])#["sif"].var[:]
end

# ╔═╡ 6288c4c6-6585-11eb-24c3-4b39125563bb
#timeStruc = ds["TIME"].attrib["units"]

# ╔═╡ 35ce9bee-62a7-11eb-0118-9528ab9f33d8
begin
	
	tropoDict = JSON.Parser.parse("""
        [{
			"basics": {
        	"lat": "lat",
        	"lon": "lon",
        	"lat_bnd": "lat_bnds",
        	"lon_bnd": "lon_bnds",
			"sif": "sif",
			"sif_err": "sif_err",
			"nir": "NIR",
        	"time": "TIME"
			},
			"time_struc": "UTC, seconds since 1970-01-01 00:00:00"
		}]
	""")

	oco2Dict = JSON.Parser.parse("""
        [{
			"basics": {
  	     	"lat": "latitude",
        	"lon": "longitude",
        	"lat_bnd": "footprint_vertex_latitude",
        	"lon_bnd": "footprint_vertex_longitude",
			"sif": "SIF_757nm",
			"sif_err": "SIF_757nm_uncert",
			"measurement_mode": "measurement_mode",
        	"time": "time",
			"nir": "continuum_radiance_757nm"
			},
			"time_struc": "seconds since 1993-1-1 0:0:0"
		}]
	""")
	
end

# ╔═╡ 254018d0-62b7-11eb-2b1d-138ad8d088b1
## TROPOMI:
#baseDict = tropoDict[1]["basics"]
## OCO-2:
baseDict = oco2Dict[1]["basics"]

# ╔═╡ 82ade7fe-62b7-11eb-3aa4-e9b8ba1f5a58
#tmp = Dataset(allFiles[1])[baseDict["lat"]].var[:]

# ╔═╡ 50f46e82-64c3-11eb-29d9-933218ab2e66
md"## Create a DataFrame from NC files"

# ╔═╡ b30830b2-71a3-11eb-2699-53e1f7bcdb43
md"
> [DataFrames.jl](https://dataframes.juliadata.org/stable/) makes it easier to work with tabular data. Its design and functionality are similar to those of pandas in Python and data.frame/data.table in R.
"

# ╔═╡ a33e67bc-650a-11eb-2ad2-43fe662d8a78
md"**Collecting data for the specified time period:**"

# ╔═╡ 7ff0299c-7265-11eb-326b-e5659d2f2f58
#filter!("measurement_mode" => ==(0), df)

# ╔═╡ 7563ebee-64f5-11eb-30ca-f5e9d7856dd7
begin
	
numberNames = function(what::String,nWhat::Int)
	numNames = Vector()
	for x in string.(Vector(1:nWhat))
		push!(numNames,what*x)
	end
	return(numNames)
end

# now return DataFrame of filtered data:
readNCdata = function(NCFiles,baseDict,mbr)
	
	dfOut = DataFrame()
	
	for aFile in NCFiles
	
		dfTmp = DataFrame()
		
		try	
		
		for aKey in keys(baseDict)
			tmp = Dataset(aFile)[baseDict[aKey]].var[:]
			if length(size(tmp))==1
				dfTmp[!,aKey] = tmp
			else
				## check which dimension is sounding and lat/lon	
				if size(tmp)[2] > size(tmp)[1] 
					extendedKeys=numberNames(aKey,size(tmp)[1])
					for iKey in 1:size(tmp)[1]
						dfTmp[!,extendedKeys[iKey]] = tmp[iKey,:]
					end				
				else
					extendedKeys=numberNames(aKey,size(tmp)[2])
					for iKey in 1:size(tmp)[2]
						dfTmp[!,extendedKeys[iKey]] = tmp[:,iKey]
					end
				end
			end
		end
			
		# filtering:
		filter!("lon" => >=(mbr[1]), dfTmp)
		filter!("lon" => <=(mbr[2]), dfTmp)
		filter!("lat" => >=(mbr[3]), dfTmp)
		filter!("lat" => <=(mbr[4]), dfTmp)		

		
		catch
		@warn "Could not read file."
		end
			
	# attach to final output data:
		if size(dfTmp)!=(0,0)
			if size(dfOut)==(0,0) 
				dfOut = dfTmp
			else
				append!(dfOut, dfTmp)
			end
		end
	
	end

	return(dfOut)

end
	
end

# ╔═╡ 00f8867e-64f3-11eb-18f0-2f31200f7628
df = readNCdata(allFiles[FileInd],baseDict,mbr)

# ╔═╡ eea1155e-64f5-11eb-3eed-f34a6000d5ea
#test numberNames function:
#numberNames("lat",4)

# ╔═╡ a3b5d552-650e-11eb-0db6-cdd12b877af8
md"**Convert time structure:**"

# ╔═╡ c21d0ec8-71a4-11eb-1aa7-b1f03c83a3b7
md" For **TROPOMI** it is unix time, so the conversion to a DateTime object is simple:"

# ╔═╡ a12da5e2-71a4-11eb-3f91-437744720d0d
tropoDict[1]["time_struc"]

# ╔═╡ eb55f0c2-6583-11eb-2572-1b96530ce77a
#df.utc = unix2datetime.(df.time)

# ╔═╡ f62cce7c-725a-11eb-3b3b-ad594257cda6
md" For **OCO-2** the time structure has to be converted in a slightly different way:"

# ╔═╡ f33daa04-725a-11eb-2831-6fc756438afe
oco2Dict[1]["time_struc"]

# ╔═╡ 724200f0-7262-11eb-1b6a-3723b249389d
begin
	timeOffset = datetime2unix(DateTime(1993,1,1, 0,0,0))
	df.utc = unix2datetime.(df.time .+ timeOffset)
end

# ╔═╡ f3550094-71a5-11eb-1bd6-8150f1ad5604
md"Add also a column with plain days for convenience:"

# ╔═╡ b99cff74-650e-11eb-3373-1d52b025e9a5
df.Date = Date.(df.utc)

# ╔═╡ 4a8de51e-65db-11eb-2e99-c70cb7193464
md"## What percentage of each footprint covers the spatial polygon?"

# ╔═╡ 26556eb4-664f-11eb-16f6-2b8c0c1a5945
## basic test for append!:
#begin
#	a = [1.]
#	#b = collect(1:10)*0.1
#	b=append!(a,(1. .+ collect(0:10) * 0.1))
#	c=append!(a,(2. .+ collect(0:10) * 0.5))
#	
#end

# ╔═╡ ac7223a8-668d-11eb-1906-a51e5b2ce8c1
md"A few required functions:"

# ╔═╡ 9555b186-664d-11eb-088a-bd3beb009eab
begin

# function to return formatted vertices:
getSdngVerts = function(df::DataFrame, i::Int64)
	# may need to edit this function later to order the vertices automatically	
	# e.g. in case LatLons are mixed up in the original data set
	vertLat = [df.lat_bnd1[i],df.lat_bnd2[i],df.lat_bnd3[i],df.lat_bnd4[i]]
	vertLon = [df.lon_bnd1[i],df.lon_bnd2[i],df.lon_bnd3[i],df.lon_bnd4[i]]
	return(DataFrame(lon=vertLon,lat=vertLat))
end
	
# To convert sounding to equally spaced spatial points we need two functions:
divLine = function(lon1,lat1,lon2,lat2,n::UInt8)
		dLat = (lat2-lat1)/(n-1)
		dLon = (lon2-lon1)/(n-1)
		
		lats = append!([lat1], (lat1 .+ collect(1:(n-1)) * dLat))
		lons = append!([lon1], (lon1 .+ collect(1:(n-1)) * dLon))
		return(lons, lats)
end
	
# Now divide each rectangle into n points:
getPoints = function(sdngVerts::DataFrame, n::UInt8)
		vertLon = sdngVerts.lon
		vertLat = sdngVerts.lat
		# Get reference points for two lines at the extremes:
	  	line1lon, line1lat = divLine(vertLon[1],vertLat[1],vertLon[2],vertLat[2],n)
	  	line2lon, line2lat = divLine(vertLon[4],vertLat[4],vertLon[3],vertLat[3],n)
	    # Now generate the inner points
		LonLat = []
		for i in 1:n
			tmpLon, tmpLat = divLine(line1lon[i], line1lat[i],  
						   			 line2lon[i], line2lat[i], n)
			for i in 1:length(tmpLon)[1]
				push!(LonLat,(tmpLon[i], tmpLat[i]))
			end
		end
		
		return(Tuple.(LonLat))
end

	
# functions to loop through soundings in DataFrame (df) to return coverage of polygon (percentage if n=10)
compCoverage = function(polygon::Array{Tuple{Float64,Float64},1}, sdngVerts::DataFrame, n::UInt8)
	testPoints = getPoints(sdngVerts, n)
	inside = [inpolygon(p, polygon; in=true, on=false, out=false) for p in testPoints]
	return(length(inside[inside .==true]))
end

getCoverage = function(df::DataFrame, polygon::Array{Tuple{Float64,Float64},1})
	[compCoverage(polygon, getSdngVerts(df,i), UInt8(10)) for i in 1:size(df)[1]]
end
	
	
end

# ╔═╡ c88c215e-6650-11eb-35fe-4dbe87cfb463
## Just to test the new functions:
# 
#getSdngVerts(df,1)
#
#divLine(Float32(0.),Float32(0.),Float32(10.), Float32(10.), UInt8(4))
#
#getPoints(getSdngVerts(df,1), UInt8(10))
#
#getCoverage(vegPolygon, getSdngVerts(df,25), UInt8(10))

# ╔═╡ a5b67f32-668f-11eb-090f-07b24fe0660d
md"**Attach coverage of all polygons to DataFrame:**"

# ╔═╡ 964d3d38-665d-11eb-1eb4-adbb37c289b6
begin
	# 
	df.VegCov = getCoverage(df,vegPolygon)
	df.RefCov = getCoverage(df,refPolygon)
end

# ╔═╡ 98684c54-6660-11eb-175f-c96664cc2f04
begin
# how does the coverage per area looks like?
#	Plots.plot(df.VegCov)
#	Plots.plot!(df.RefCov)
end

# ╔═╡ 0ff6e15e-6661-11eb-3e48-3b73b977b0fb
md"## Let's take a first look at the data:"

# ╔═╡ 0cda4bd4-6691-11eb-20c7-9b5097e1e7e9
begin
	minCovVeg_slider= @bind minCovVeg Slider(0:10:100; default=50)
	minCovRef_slider= @bind minCovRef Slider(0:10:100; default=50)
	md"""
	**Vegetation:** $(minCovVeg_slider)
	
	**Reference:**  $(minCovRef_slider)
	"""
end

# ╔═╡ 388e8398-6692-11eb-16c8-2bcccd929b61
md"**Soundings in time series have a coverage of at least $(minCovVeg)% (Vegetation) and $(minCovRef)% (Reference)**"

# ╔═╡ 9ba7d79a-6692-11eb-30e8-8d548156f71c
begin
	indVeg  = findall(df.VegCov .>= minCovVeg)
	ntotVeg = length(findall(df.VegCov .> 0))#size(df)[1] 

	Plots.scatter(df.utc[indVeg],df.sif[indVeg], markeralpha=0.5, label= "Vegetation (n="*string(length(indVeg))*"/"*string(ntotVeg)*")")

	indRef = findall(df.RefCov .>= minCovRef)
	ntotRef = length(findall(df.RefCov .> 0))#size(df)[1] 
	Plots.scatter!(df.utc[indRef],df.sif[indRef], markeralpha=0.5, label= "Reference (n="*string(length(indRef))*"/"*string(ntotRef)*")")
	
end

# ╔═╡ 46c0b6a0-6824-11eb-2407-a37dca8ef4bf
md"## Temporal averaging"

# ╔═╡ 2140d9a2-682c-11eb-0be0-8d4cba6681bb
begin
	time_slider= @bind n_days Slider(1:1:31; default=7)
	md"""
	**Averaging period (days):** $(time_slider)
	"""
end

# ╔═╡ 6988f26c-682c-11eb-3261-ed8a7c839fb4
md"**$(n_days) Day(s)** "

# ╔═╡ e28b9fc0-6a57-11eb-183f-6f47f1c39d6a
begin
#df.Date = Date.(df.utc)
	
## function returns temporal average (mean, std, n)	
tempAv = function(df::DataFrame, nDays, what::String)
	
	# rename the column of interest
	subDF = select(df, what => "what", "Date" => "Date")	
	
	# generate breaks:	
	tBreaks = collect(minimum(subDF.Date):Day(nDays):maximum(subDF.Date)+Day(1))
	
	## initialize avTime column:
	subDF.avTime = Array{DateTime, 1}(undef,  size(subDF)[1])
	
	## populate avTime:
	for i in 1:(length(tBreaks)-1)
			subDF.avTime[findall((subDF.Date .>= tBreaks[i]) .& (subDF.Date .< tBreaks[i+1]))] .= DateTime(tBreaks[i]) .+ (DateTime(tBreaks[i+1]) - DateTime(tBreaks[i])) ./ 2
	end
	## clean-up incomplete time steps:
	#filter!(row -> row[:avTime]!= DateTime("0000-12-31T00:00:00"), subDF)
	#filter!(row -> row[:avTime] > minimum(tBreaks), subDF)
	#filter!(row -> row[:avTime] < maximum(tBreaks), subDF)
	filter!("avTime" => >(minimum(tBreaks)), subDF)
	filter!("avTime" => <(maximum(tBreaks)), subDF)

	## generate output	
	outDF = combine(groupby(subDF, :avTime), :what => mean)
	rename!(outDF, :what_mean => :mean)
	outDF.sd = combine(groupby(subDF, :avTime), :what => std)[:,2]
	outDF.n   = combine(groupby(subDF, :avTime), :what => length)[:,2]
	#sort!(outDF, :avTime)
	return(outDF)	
end
	
end

# ╔═╡ 5623b812-6825-11eb-091d-3b676411207c
begin
av = tempAv(df[indVeg,:], n_days, "sif")	
scatter(df.utc[indVeg],df.sif[indVeg], markeralpha=0.2, c=:gray, label="single soundings" )

plot!(av.avTime,av.mean; ribbon=av.sd, linewidth=4, label= "mean and sd ("*string(n_days)*" days)")
end

# ╔═╡ d51c868a-760d-11eb-3420-9b51e5d9299a
av

# ╔═╡ e53f342c-669a-11eb-3afb-fbeede548218
md"## Computing spatial composites"

# ╔═╡ 59c35358-682d-11eb-0fd2-77744637cdbe
md"""**Resolution:**
		$(@bind resolution Select(["0.2" => "0.2°","0.1" => "0.1°","0.05" => "0.05°", "0.01" => "0.01°"]))
		"""

# ╔═╡ 2e2828e8-67ee-11eb-1b37-3ffc4a162aeb
begin
	
initGrid = function(resolution::Float64, mbr::Array{Float64,1}) ##mbr format: lonmin,lonmax,latmin,latmax
		
	lons = collect(mbr[1]+resolution/2.:resolution:mbr[2]-resolution/2.)
	lats = collect(mbr[3]+resolution/2.:resolution:mbr[4]-resolution/2.)
	
    #lons = collect(mbr[1]:resolution:mbr[2])
	#lats = collect(mbr[3]:resolution:mbr[4])
	
	iLons = collect(1:length(lons))
	iLats = collect(1:length(lats))
	
	## to map lon/lat to iLon/iLat
	iLon_ = Spline1D(lons,iLons)
	iLat_ = Spline1D(lats,iLats)

	grid = Array{Union{Missing,Float32}}(undef, length(lats), length(lons))
	grid[:,:] .= 0.
	
	return(lons, lats, iLon_, iLat_, grid)	
end

## generate polygon from corner coordinates:	
getPoly = function(xs,ys)
##first and last Tuple neeed to be identical to close the polygon
		polygon = []
		for i in 1:length(xs)
			push!(polygon, (xs[i], ys[i]))
		end
		push!(polygon, (xs[1], ys[1]))
	
		return(Tuple.(polygon))
end

	
oversample = function(df::DataFrame, resolution::Float64, mbr::Array{Float64,1}, what::String)
	
	#initialize grid:	
	lons, lats, iLon_,iLat_,grid = initGrid(resolution, mbr)	
	nGrid = copy(grid)
	
	## loop through soundings:	
	for iSdng in 1:size(df)[1]
		
		sdngVert = getSdngVerts(df,iSdng)
		
		# convert from coordinates to XY
		xs = iLon_(sdngVert.lon) 
		ys = iLat_(sdngVert.lat)	
		
		polygon = getPoly(xs,ys)
		
		# loop through all gridboxes in that area: 
		xmin = Int64(floor(minimum(xs)))
		if xmin==0 
			xmin += 1
		end
		xmax = Int64(ceil(maximum(xs)))
		if xmax > length(lons)
			xmax -= 1
		end
		ymin = Int64(floor(minimum(ys)))
		if ymin==0 
			ymin += 1
		end
		ymax = Int64(ceil(maximum(ys)))
		if ymax > length(lats)
			ymax -= 1
		end
			
		# compute coverage only if grid box vertices are not covered by sounding
		for x in xmin:xmax
			for y in ymin:ymax
				
				testPoints = [(x,y),(x,y+1),(x+1,y+1),(x+1,y)] # clockwise
				
				inside = [inpolygon(p, polygon; in=true, on=true, out=false) for p in testPoints]	
				## easy case if all corners are within sounding:
				if length(inside[inside .==true])==4
					nGrid[y,x] += 1.
					grid[y,x]  += df[iSdng,what]
				else	
				## compute coverage if gridcell is only partially covered:
					inflatedPoints = getPoints(DataFrame(lon=Float32.([x, x, x+1, x+1]), lat=Float32.([y, y+1, y+1, y])), UInt8(10))
					cov = [inpolygon(p, polygon; in=true, on=true, out=false) for p in inflatedPoints]
					weight = Float32(length(cov[cov .==true])/100.)
					nGrid[y,x] += weight
					grid[y,x]  += df[iSdng,what] * weight
					
				end
			end
		end
	end
	nGrid[nGrid .==0] .= missing
	avGrid = grid ./ nGrid

	return(avGrid, lons, lats)
end
	
	
end

# ╔═╡ 4a72d39c-65da-11eb-07c2-2dd1c2bcac1e
begin
	### use the first sounding with coverage above 10%:
	id = findfirst((df.VegCov .> 20.).&(df.VegCov .< 30.))
	
	Plots.plot(vegPolygon, legend=false, fill = (34,0.5,:gray))
	Sdng = getSdngVerts(df,id)
	testPoints = getPoints(Sdng , UInt8(10))
	Plots.scatter!(testPoints)
	
	touchPoly = getPoly(Sdng.lon,Sdng.lat)
	Plots.plot!(touchPoly, linewidth=4)
	title!("Coverage: "*string(df.VegCov[id])*"%")
end

# ╔═╡ 620677cc-6817-11eb-0205-e9daaf8aaa37
begin
let
	## use just one month of data for faster gridding:
	idxTime = findall(month.(df.Date) .== 5)
	## or two month?:
	#in_month = in([4,5,6])
	#idxTime = findall(in_month.(month.(df.Date)).==true)
	## Or half a month?:
	#idxTime = findall((month.(df.Date) .== 5) .& (day.(df.Date) .<= 15))

	out, lons, lats = oversample(df[:,:],tryparse(Float64,resolution), 										 mbr,"sif")
	
		
	Plots.heatmap(lons, lats, out, color= :viridis)#, clim=(-0.2,0.9))
	#fin = oversample(df, resolution, mbr,"sif")
	#Plots.plot(fin, color= :viridis)
	Plots.plot!(vegPolygon, legend=false, linewidth=3)
	Plots.plot!(refPolygon, legend=false, linewidth=3, linecolor=:blue)
	title!("n="*string(length(idxTime)))
end
end

# ╔═╡ f1f9738c-669b-11eb-3a6f-011326f23714
md"## Uncertainty/measurement noise:"

# ╔═╡ d1cec2be-71a7-11eb-1d36-e9c097be9b7f
md"
> **Negative SIF values are plausible due to the measurement uncertainty. Filtering negative SIF values will result in a positive bias!**
"

# ╔═╡ a9471dee-6b0a-11eb-0dde-19829d25aab5
begin

	#rnorm(n, mu, sig) = rand(Normal(mu, sig), n)
	
	### what can be expected (mean error estimate, uncertainty prediction):
	meanErr = mean(df[indRef,"sif_err"])
	
	### what is actually observed:
	fitErr = fit(Normal, df[indRef,"sif"])

	## plot:
	histogram(df[indRef,"sif"], bins = 100, normalize = :pdf, label="Observed σ = "*string(round(fitErr.σ; digits=2))*", μ = "*string(round(fitErr.μ; digits=2)))
	plot!(Normal(0.,meanErr), linewidth=4, label="Predicted σ = "*string(round(meanErr; digits=2))*", μ = 0")	
	title!("n="*string(length(indRef)))

end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
Dierckx = "39dd38d3-220a-591b-8e3c-4c3a8c710a94"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
GeoInterface = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
GeoJSON = "61d90e0f-e114-555e-ac52-39dfb47a3ef9"
Glob = "c27321d9-0574-5035-807b-f59d2c89b15c"
HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"
JSON = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
NCDatasets = "85f8d34a-cbdd-5861-8df4-14fed0d494ab"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PolygonOps = "647866c9-e3ac-4575-94e7-e3d426903924"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsPlots = "f3b207a7-027a-5e70-b257-86293d7955fd"

[compat]
DataFrames = "~1.2.2"
Dierckx = "~0.5.1"
Distributions = "~0.25.11"
GeoInterface = "~0.5.5"
GeoJSON = "~0.5.1"
Glob = "~1.3.0"
HTTP = "~0.9.13"
JSON = "~0.21.2"
NCDatasets = "~0.11.7"
PlutoUI = "~0.7.9"
PolygonOps = "~0.1.1"
StatsPlots = "~0.14.26"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "485ee0867925449198280d4af84bdb46a2a404d0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.0.1"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra"]
git-tree-sha1 = "2ff92b71ba1747c5fdd541f8fc87736d82f40ec9"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.4.0"

[[Arpack_jll]]
deps = ["Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "e214a9b9bd1b4e1b4f15b22c0994862b66af7ff7"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.0+3"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "a4d07a1c313392a77042855df46c5f534076fab9"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.0"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c3598e525718abcc440f69cc6d5f60dda0a1b61e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.6+5"

[[CFTime]]
deps = ["Dates", "Printf"]
git-tree-sha1 = "bca6cb6ee746e6485ca4535f6cc29cf3579a0f20"
uuid = "179af706-886a-5703-950a-314cd64e0468"
version = "0.1.1"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "e2f47f6d8337369411569fd45ae5753ca10394c6"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.0+6"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "bdc0937269321858ab2a4f288486cb258b9a0af7"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.3.0"

[[Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "75479b7df4167267d75294d14b58244695beb2ac"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.14.2"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "9995eb3977fbf67b86d0a0a0508e83017ded03f2"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.14.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "727e463cfebd0c7b999bbf3e9e7e16f254b94193"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.34.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[Crayons]]
git-tree-sha1 = "3f71217b538d7aaee0b69ab47d9b7724ca8afa0d"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.0.4"

[[DataAPI]]
git-tree-sha1 = "ee400abb2298bd13bfc3df1c412ed228061a2385"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.7.0"

[[DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d785f42445b63fc86caa08bb9a9351008be9b765"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.2.2"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[DataValues]]
deps = ["DataValueInterfaces", "Dates"]
git-tree-sha1 = "d88a19299eba280a6d062e135a43f00323ae70bf"
uuid = "e7dc6d0d-1eca-5fa6-8ad6-5aecde8b7ea5"
version = "0.4.13"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Dierckx]]
deps = ["Dierckx_jll"]
git-tree-sha1 = "5fefbe52e9a6e55b8f87cb89352d469bd3a3a090"
uuid = "39dd38d3-220a-591b-8e3c-4c3a8c710a94"
version = "0.5.1"

[[Dierckx_jll]]
deps = ["CompilerSupportLibraries_jll", "Libdl", "Pkg"]
git-tree-sha1 = "a580560f526f6fc6973e8bad2b036514a4e3b013"
uuid = "cd4c43a9-7502-52ba-aa6d-59fb2a88580b"
version = "0.0.1+0"

[[Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "abe4ad222b26af3337262b8afb28fab8d215e9f8"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.3"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns"]
git-tree-sha1 = "3889f646423ce91dd1055a76317e9a1d3a23fff1"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.11"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "a32185f5428d3986f47c2ab78b1f216d5e6cc96f"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.5"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "92d8f9f208637e8d2d28c664051a00569c01493d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.1.5+1"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "LibVPX_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "3cc57ad0a213808473eafef4845a74766242e05f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.3.1+4"

[[FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "f985af3b9f4e278b1d24434cbb546d6092fca661"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.3"

[[FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3676abafff7e4ff07bbd2c42b3d8201f31653dcc"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.9+8"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "7c365bdef6380b29cfc5caaf99688cd7489f9b87"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.2"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "35895cf184ceaab11fd778b4590144034a167a2f"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.1+14"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "cbd58c9deb1d304f5a245a0b7eb841a2560cfec6"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.1+5"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "dba1e8614e98949abfa60480b13653813d8f0157"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+0"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "182da592436e287758ded5be6e32c406de3a2e47"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.58.1"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d59e8320c2747553788e4fc42231489cc602fa50"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.58.1+0"

[[GeoInterface]]
deps = ["RecipesBase"]
git-tree-sha1 = "38a649e6a52d1bea9844b382343630ac754c931c"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "0.5.5"

[[GeoJSON]]
deps = ["GeoInterface", "JSON3"]
git-tree-sha1 = "4764da92d333658552b2bedc9f6b379f017c727b"
uuid = "61d90e0f-e114-555e-ac52-39dfb47a3ef9"
version = "0.5.1"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "7bf67e9a481712b3dbe9cb3dac852dc4b1162e02"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+0"

[[Glob]]
git-tree-sha1 = "4df9f7e06108728ebf00a0a11edee4b29a482bb2"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HDF5_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "OpenSSL_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "fd83fa0bde42e01952757f01149dd968c06c4dba"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.12.0+1"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "44e3b40da000eab4ccb1aecdc4801c040026aeb5"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.13"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "61aa005707ea2cebf47c8d780da8dc9bc4e0c512"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.4"

[[InvertedIndices]]
deps = ["Test"]
git-tree-sha1 = "15732c475062348b0165684ffe28e85ea8396afc"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.0.0"

[[IrrationalConstants]]
git-tree-sha1 = "f76424439413893a832026ca355fe273e93bce94"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.0"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JSON3]]
deps = ["Dates", "Mmap", "Parsers", "StructTypes", "UUIDs"]
git-tree-sha1 = "b3e5984da3c6c95bcf6931760387ff2e64f508f3"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.9.1"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "591e8dc09ad18386189610acafb970032c519707"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.3"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a4b12a1bd2ebade87891ab7e36fdbce582301a92"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.6"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[LibVPX_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "12ee7e23fa4d18361e7c2cde8f8337d4c3101bc7"
uuid = "dd192d2f-8180-539f-9fb4-cc70b1dcf69a"
version = "1.10.0+0"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "761a393aeccd6aa92ec3515e428c26bf99575b3b"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+0"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "3d682c07e6dd250ed082f883dc88aee7996bf2cc"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.0"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "c253236b0ed414624b083e6b72bfe891fbd2c7af"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2021.1.1+1"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "0fb723cd8c45858c22169b2e42269e53271a6df7"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.7"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "2ca267b08821e86c5ef4376cffed98a46c2cb205"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.1"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MultivariateStats]]
deps = ["Arpack", "LinearAlgebra", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "8d958ff1854b166003238fe191ec34b9d592860a"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.8.0"

[[NCDatasets]]
deps = ["CFTime", "DataStructures", "Dates", "NetCDF_jll", "Printf"]
git-tree-sha1 = "5da406d9624f25909a6f556bd8d5c1deaa189ee6"
uuid = "85f8d34a-cbdd-5861-8df4-14fed0d494ab"
version = "0.11.7"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "16baacfdc8758bc374882566c9187e785e85c2f0"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.9"

[[NetCDF_jll]]
deps = ["Artifacts", "HDF5_jll", "JLLWrappers", "LibCURL_jll", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Pkg", "Zlib_jll", "nghttp2_jll"]
git-tree-sha1 = "0cf4d1bf2ef45156aed85c9ac5f8c7e697d9288c"
uuid = "7243133f-43d8-5620-bbf4-c2c921802cf3"
version = "400.702.400+0"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Observables]]
git-tree-sha1 = "fe29afdef3d0c4a8286128d4e45cc50621b1e43d"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.4.0"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "c0f4a4836e5f3e0763243b8324200af6d0e0f90c"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.5"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7937eda4681660b4d6aeeecc2f7e1c81c8ee4e2f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+0"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "4dd403333bcf0909341cfe57ec115152f937d7d8"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.1"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "438d35d2d95ae2c5e8780b330592b6de8494e779"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.3"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "501c20a63a34ac1d015d5304da0e645f42d91c9f"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.11"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "8365fa7758e2e8e4443ce866d6106d8ecbb4474e"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.20.1"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[PolygonOps]]
git-tree-sha1 = "c031d2332c9a8e1c90eca239385815dc271abb22"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.1"

[[PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "cde4ce9d6f33219465b55162811d8de8139c0414"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.2.1"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "0d1245a357cc61c8cd61934c07447aa569ff22e6"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.1.0"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "12fbe86da16df6679be7521dfb39fbc861e1dc7b"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.1"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Ratios]]
deps = ["Requires"]
git-tree-sha1 = "7dff99fbc740e2f8228c6878e2aad6d7c2678098"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.1"

[[RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "2a7a2469ed5d94a98dea0e85c46fa653d76be0cd"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.3.4"

[[Reexport]]
git-tree-sha1 = "adcd36e8ba9665c88eb0bd156d4e2a19f9b0d889"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.0"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "54f37736d8934a12a200edea2f9206b03bdf3159"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.7"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "LogExpFunctions", "OpenSpecFun_jll"]
git-tree-sha1 = "a322a9493e49c5f3a10b50df3aedaf1cdb3244b7"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.6.1"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3240808c6d463ac46f1c1cd7638375cd22abbccb"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.12"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "fed1ec1e65749c4d96fc20dd13bea72b55457e62"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.9"

[[StatsFuns]]
deps = ["IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "20d1bb720b9b27636280f751746ba4abb465f19d"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.9"

[[StatsPlots]]
deps = ["Clustering", "DataStructures", "DataValues", "Distributions", "Interpolations", "KernelDensity", "LinearAlgebra", "MultivariateStats", "Observables", "Plots", "RecipesBase", "RecipesPipeline", "Reexport", "StatsBase", "TableOperations", "Tables", "Widgets"]
git-tree-sha1 = "e7d1e79232310bd654c7cef46465c537562af4fe"
uuid = "f3b207a7-027a-5e70-b257-86293d7955fd"
version = "0.14.26"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "000e168f5cc9aded17b6999a560b7c11dda69095"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.0"

[[StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "e36adc471280e8b346ea24c5c87ba0571204be7a"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.7.2"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableOperations]]
deps = ["SentinelArrays", "Tables", "Test"]
git-tree-sha1 = "019acfd5a4a6c5f0f38de69f2ff7ed527f1881da"
uuid = "ab02a1b2-a7df-11e8-156e-fb1833f50b87"
version = "1.1.0"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "d0c690d37c73aeb5ca063056283fde5585a41710"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.5.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll"]
git-tree-sha1 = "2839f1c1296940218e35df0bbb220f2a79686670"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.18.0+4"

[[Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "eae2fbbc34a79ffd57fb4c972b08ce50b8f6a00d"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.3"

[[WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "59e2ad8fd1591ea019a5259bd012d7aee15f995c"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.3"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "acc685bcf777b2202a904cdcb49ad34c2fa1880c"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.14.0+4"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7a5780a0d9c6864184b3a2eeeb833a0c871f00ab"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "0.1.6+4"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d713c1ce4deac133e3334ee12f4adff07f81778f"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2020.7.14+2"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "487da2f8f2f0c8ee0e83f39d13037d6bbf0a45ab"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.0.0+3"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─6659c6d4-7182-11eb-3b65-d939a40b4325
# ╟─1c2705ee-61a4-11eb-33f9-835cc5caa8a8
# ╟─b7921386-7186-11eb-20d2-d973f7940a35
# ╟─e7799730-7187-11eb-02a8-e7f76ed6081c
# ╟─b434a196-61a4-11eb-19fe-976e23dec488
# ╠═c25b5f6c-61a4-11eb-2c34-cbf860b5d314
# ╟─a4854224-7184-11eb-291a-e3667b3c33e2
# ╠═5250d5d2-718f-11eb-03ab-5d29ef4ce6ad
# ╟─34be05d8-7185-11eb-3bd3-d1e24885ed75
# ╟─8cc541f6-7185-11eb-3034-675a8bf6ff7c
# ╟─2be873b4-718d-11eb-1a33-a742dd4610ff
# ╟─bb5cb8a2-718d-11eb-3694-cbe3adf78a8d
# ╟─a33e9ac2-718c-11eb-18d4-b90c627f52ec
# ╟─8f48bb84-718e-11eb-333e-8f1f2df863d9
# ╟─38f5b8de-7190-11eb-0437-497aa85610f0
# ╟─2df3dee4-61a4-11eb-2003-5b56cfb49d01
# ╠═5b5b38b4-61a4-11eb-284b-e329267e0fb2
# ╟─4f00b9e0-61a4-11eb-0d33-f3ec82cd5396
# ╠═2f0d9a3a-61a5-11eb-2fc1-650e465800c8
# ╟─62fa453c-61a5-11eb-0b39-015a9ff4155c
# ╠═a18e5c16-61c3-11eb-37cf-edd1266f9dde
# ╠═179972c6-7257-11eb-3a3d-45c672612fb3
# ╟─bc20f27e-61d1-11eb-2497-179d53ab2360
# ╠═dfaee354-61d1-11eb-2a88-8bb2002c2136
# ╟─a8dc819e-7188-11eb-2c1b-8998b48d99f5
# ╟─cfd6f918-7189-11eb-33f4-c5133381580b
# ╟─c475fc80-718a-11eb-2c79-3f45c09aaca4
# ╟─2680b0a4-718c-11eb-2bad-137068539830
# ╟─471fe14a-718c-11eb-23a2-a3c5ed97e288
# ╟─a748c364-71a2-11eb-21fb-2ddc00692957
# ╟─a8564f6a-71a2-11eb-1797-67abf06ac8c1
# ╟─d753c838-71a2-11eb-151d-6394549ab869
# ╟─d7f995ce-71a2-11eb-3e35-b113a803a9fb
# ╟─91b73620-7193-11eb-070f-573411f9e611
# ╠═df1138ee-61d5-11eb-34b6-67ce049efb74
# ╠═14e5f60c-6687-11eb-37f6-fbf0d9061530
# ╠═ad814292-65d9-11eb-1390-35f78905835c
# ╠═e083f3fc-65da-11eb-290c-8df8464a4ff1
# ╟─53d0a158-62a6-11eb-3f8b-af81223169ec
# ╟─213784d6-7198-11eb-16a2-254812fbed1d
# ╟─b88cdc64-6509-11eb-24c7-fd78be2de33c
# ╠═80cc3a5c-629f-11eb-3807-f51439287448
# ╠═b1b2bb48-64fd-11eb-1ce8-274175a2d8ff
# ╟─b3fa469c-61d4-11eb-17ca-31a9319d9560
# ╠═50938100-62af-11eb-0027-2bc902261db2
# ╟─6288c4c6-6585-11eb-24c3-4b39125563bb
# ╠═35ce9bee-62a7-11eb-0118-9528ab9f33d8
# ╠═254018d0-62b7-11eb-2b1d-138ad8d088b1
# ╟─82ade7fe-62b7-11eb-3aa4-e9b8ba1f5a58
# ╟─50f46e82-64c3-11eb-29d9-933218ab2e66
# ╟─b30830b2-71a3-11eb-2699-53e1f7bcdb43
# ╟─a33e67bc-650a-11eb-2ad2-43fe662d8a78
# ╠═00f8867e-64f3-11eb-18f0-2f31200f7628
# ╠═7ff0299c-7265-11eb-326b-e5659d2f2f58
# ╠═7563ebee-64f5-11eb-30ca-f5e9d7856dd7
# ╟─eea1155e-64f5-11eb-3eed-f34a6000d5ea
# ╟─a3b5d552-650e-11eb-0db6-cdd12b877af8
# ╟─c21d0ec8-71a4-11eb-1aa7-b1f03c83a3b7
# ╠═a12da5e2-71a4-11eb-3f91-437744720d0d
# ╠═eb55f0c2-6583-11eb-2572-1b96530ce77a
# ╟─f62cce7c-725a-11eb-3b3b-ad594257cda6
# ╠═f33daa04-725a-11eb-2831-6fc756438afe
# ╠═724200f0-7262-11eb-1b6a-3723b249389d
# ╟─f3550094-71a5-11eb-1bd6-8150f1ad5604
# ╠═b99cff74-650e-11eb-3373-1d52b025e9a5
# ╟─4a8de51e-65db-11eb-2e99-c70cb7193464
# ╟─26556eb4-664f-11eb-16f6-2b8c0c1a5945
# ╠═4a72d39c-65da-11eb-07c2-2dd1c2bcac1e
# ╟─ac7223a8-668d-11eb-1906-a51e5b2ce8c1
# ╠═9555b186-664d-11eb-088a-bd3beb009eab
# ╟─c88c215e-6650-11eb-35fe-4dbe87cfb463
# ╟─a5b67f32-668f-11eb-090f-07b24fe0660d
# ╠═964d3d38-665d-11eb-1eb4-adbb37c289b6
# ╟─98684c54-6660-11eb-175f-c96664cc2f04
# ╟─0ff6e15e-6661-11eb-3e48-3b73b977b0fb
# ╠═0cda4bd4-6691-11eb-20c7-9b5097e1e7e9
# ╠═388e8398-6692-11eb-16c8-2bcccd929b61
# ╠═9ba7d79a-6692-11eb-30e8-8d548156f71c
# ╟─46c0b6a0-6824-11eb-2407-a37dca8ef4bf
# ╟─6988f26c-682c-11eb-3261-ed8a7c839fb4
# ╟─2140d9a2-682c-11eb-0be0-8d4cba6681bb
# ╠═5623b812-6825-11eb-091d-3b676411207c
# ╠═d51c868a-760d-11eb-3420-9b51e5d9299a
# ╠═e28b9fc0-6a57-11eb-183f-6f47f1c39d6a
# ╟─e53f342c-669a-11eb-3afb-fbeede548218
# ╟─59c35358-682d-11eb-0fd2-77744637cdbe
# ╠═620677cc-6817-11eb-0205-e9daaf8aaa37
# ╠═2e2828e8-67ee-11eb-1b37-3ffc4a162aeb
# ╟─f1f9738c-669b-11eb-3a6f-011326f23714
# ╟─d1cec2be-71a7-11eb-1d36-e9c097be9b7f
# ╠═a9471dee-6b0a-11eb-0dde-19829d25aab5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
