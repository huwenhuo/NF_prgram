#!/usr/bin/env python3
import os
import re

def split_nextflow_modules(input_file, output_dir="modules"):
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    with open(input_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Regex pattern to match an entire process block
    # Looks for process PROCESS_NAME { ... } spanning multiple lines
    process_pattern = re.compile(
        r'(process\s+([A-Z0-9_]+)\s*\{.*?\n\})', 
        re.DOTALL
    )
    
    matches = process_pattern.findall(content)
    
    if not matches:
        print("No processes found! Ensure your process blocks end with '}' on a new line.")
        return

    print(f"Found {len(matches)} processes in '{input_file}'. Splitting...\n")
    
    include_statements = []

    for full_process_code, process_name in matches:
        # Convert PROCESS_NAME to lowercase for file naming (e.g., DESEQ2_IRFINDER -> deseq2_irfinder.nf)
        file_name = f"{process_name.lower()}.nf"
        file_path = os.path.join(output_dir, file_name)
        
        # Write individual process file
        with open(file_path, 'w', encoding='utf-8') as out_f:
            out_f.write(full_process_code.strip() + "\n")
            
        print(f" Saved: {file_path}")
        
        # Generate corresponding DSL2 include statement for your workflow
        rel_path = f"./{output_dir}/{process_name.lower()}"
        include_statements.append(f"include {{ {process_name} }} from '{rel_path}'")

    print("\n" + "="*50)
    print("FINISHED! Add these includes to your main.nf or workflow file:")
    print("="*50 + "\n")
    for stmt in include_statements:
        print(stmt)

if __name__ == "__main__":
    split_nextflow_modules("modules.nf")
