#!/usr/bin/env python3
"""
Merge multiple CSV files from results/source_data/ into a single Excel file.
Each CSV file becomes a separate sheet in the Excel file, with the filename as the sheet name.
"""

import os
import pandas as pd
from pathlib import Path
import sys

def merge_csv_to_excel():
    """
    Merge all CSV files in results/source_data/ into a single Excel file.
    """
    # Get the project root directory (parent of src)
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    
    # Define paths
    source_data_dir = project_root / "results" / "source_data"
    output_file = project_root / "results" / "source_data_merged.xlsx"
    
    # Check if source directory exists
    if not source_data_dir.exists():
        print(f"错误: 源数据目录不存在: {source_data_dir}")
        sys.exit(1)
    
    # Find all CSV files
    csv_files = list(source_data_dir.glob("*.csv"))
    
    if not csv_files:
        print(f"警告: 在 {source_data_dir} 中没有找到 CSV 文件")
        sys.exit(1)
    
    print(f"找到 {len(csv_files)} 个 CSV 文件")
    
    # Create Excel writer object
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        processed_count = 0
        failed_files = []
        
        for csv_file in sorted(csv_files):
            try:
                # Get sheet name from filename (remove .csv extension)
                sheet_name = csv_file.stem
                
                # Excel sheet names have a 31 character limit
                if len(sheet_name) > 31:
                    sheet_name = sheet_name[:31]
                
                # Read CSV file
                df = pd.read_csv(csv_file)
                
                # Write to Excel sheet
                df.to_excel(writer, sheet_name=sheet_name, index=False)
                
                print(f"已处理: {csv_file.name} -> sheet: {sheet_name}")
                processed_count += 1
                
            except Exception as e:
                print(f"错误处理文件 {csv_file.name}: {str(e)}")
                failed_files.append(csv_file.name)
                continue
    
    # Print summary
    print(f"\n合并完成!")
    print(f"成功处理: {processed_count} 个文件")
    if failed_files:
        print(f"失败的文件: {len(failed_files)} 个")
        for file in failed_files:
            print(f"  - {file}")
    print(f"输出文件: {output_file}")
    
    return output_file

if __name__ == "__main__":
    merge_csv_to_excel()
