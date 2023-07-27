# Create header for C file
file(WRITE ${OUTPUT_C} "/* Autogenerated by bin2c */\n\n")
file(APPEND ${OUTPUT_C} "#include <stdint.h>\n\n")

# Create header of H file
file(WRITE ${OUTPUT_H} "/* Autogenerated by bin2c */\n\n")
file(APPEND ${OUTPUT_H} "#pragma once\n")
file(APPEND ${OUTPUT_H} "#include <stdint.h>\n\n")

string(REPLACE "," ";" INPUT_LIST ${INPUT_FILES})


# Iterate through binary files files
foreach(bin ${INPUT_LIST})
  # Get short filename
  string(REGEX MATCH "([^/]+)$" filename ${bin})
  # Replace filename spaces & extension separator for C compatibility
  string(REGEX REPLACE "\\.| |-" "_" filename ${filename})
  # Convert to lower case
  string(TOLOWER ${filename} filename)
  # Read hex data from file
  file(READ ${bin} filedata HEX)
  # Convert hex data for C compatibility
  string(REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1," filedata ${filedata})
  # Append data to c file
  file(APPEND ${OUTPUT_C} "uint8_t ${filename}[] = {${filedata}};\n\nuint32_t ${filename}_size = sizeof(${filename});\n\n")
  # Append extern definitions to h file
  file(APPEND ${OUTPUT_H} "extern uint8_t ${filename}[];\n\nextern uint32_t ${filename}_size;\n\n")
endforeach()