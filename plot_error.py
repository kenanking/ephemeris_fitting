"""
File: plot_error.py
File Created: Monday, 21st November 2022 11:49:35 am
Author: Yan Tang (360383464@qq.com)
-----
Last Modified: Monday, 21st November 2022 9:47:03 pm
Modified By: Yan Tang (360383464@qq.com>)
-----
Copyright 2022 - 2022 Yan Tang
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt

from argparse import ArgumentParser

# 显示中文
plt.rcParams["font.sans-serif"] = ["SimHei"]
plt.rcParams["axes.unicode_minus"] = False


def parse_args():
    parser = ArgumentParser(description="plot error")
    parser.add_argument("--result", "-r", type=str, default="result", help="result dir")
    parser.add_argument("--output", "-o", type=str, default="output", help="output dir")
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    result_dir = args.result
    output_dir = args.output

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 确认result_file是个文件
    if not os.path.isdir(result_dir):
        print("result_dir is not a dir")
        return

    files = os.listdir(result_dir)
    files = [os.path.join(result_dir, file) for file in files]

    prn_list = []
    sigma_0_list = []
    rms_x_list = []
    rms_y_list = []
    rms_z_list = []

    for f in files:
        # 读取result_file
        with open(f, "r") as f:
            data = json.load(f)

        prn_list.append(data["prn"])
        sigma_0_list.append(data["sigma_0"] * 1000)
        rms_x_list.append(data["interpolation"]["rms_x"] * 1000)
        rms_y_list.append(data["interpolation"]["rms_y"] * 1000)
        rms_z_list.append(data["interpolation"]["rms_z"] * 1000)

        interp_error = data["interpolation"]
        extrap_left_error = data["left_extrapolation"]
        extrap_right_error = data["right_extrapolation"]

        # 绘制内插误差
        plt.figure()
        plt.plot(
            interp_error["toe"],
            np.array(interp_error["x"]) * 1000,
            c="#3498DB",
            label="x error",
        )
        plt.plot(
            interp_error["toe"],
            np.array(interp_error["y"]) * 1000,
            c="#1ABC9C",
            label="y error",
        )
        plt.plot(
            interp_error["toe"],
            np.array(interp_error["z"]) * 1000,
            c="#F39C12",
            label="z error",
        )

        plt.plot(
            extrap_left_error["toe"],
            np.array(extrap_left_error["x"]) * 1000,
            c="#3498DB",
            linestyle="--",
        )
        plt.plot(
            extrap_left_error["toe"],
            np.array(extrap_left_error["y"]) * 1000,
            c="#1ABC9C",
            linestyle="--",
        )
        plt.plot(
            extrap_left_error["toe"],
            np.array(extrap_left_error["z"]) * 1000,
            c="#F39C12",
            linestyle="--",
        )

        plt.plot(
            extrap_right_error["toe"],
            np.array(extrap_right_error["x"]) * 1000,
            c="#3498DB",
            linestyle="-.",
        )
        plt.plot(
            extrap_right_error["toe"],
            np.array(extrap_right_error["y"]) * 1000,
            c="#1ABC9C",
            linestyle="-.",
        )
        plt.plot(
            extrap_right_error["toe"],
            np.array(extrap_right_error["z"]) * 1000,
            c="#F39C12",
            linestyle="-.",
        )

        # 绘制最大最小内插误差限
        max_interp_error = (
            max(max(interp_error["x"]), max(interp_error["y"]), max(interp_error["z"]))
            * 1000
        )
        plt.hlines(
            max_interp_error,
            interp_error["toe"][0],
            interp_error["toe"][-1],
            colors="#E74C3C",
            linestyles="dashed",
        )
        plt.text(
            interp_error["toe"][-1] - 1000,
            max_interp_error + 2,
            f"{max_interp_error:.2f}",
            color="#E74C3C",
            fontsize=12,
        )

        min_interp_error = (
            min(min(interp_error["x"]), min(interp_error["y"]), min(interp_error["z"]))
            * 1000
        )
        plt.hlines(
            min_interp_error,
            interp_error["toe"][0],
            interp_error["toe"][-1],
            colors="#E74C3C",
            linestyles="dashed",
        )
        plt.text(
            interp_error["toe"][-1] - 1000,
            min_interp_error - 8,
            f"{min_interp_error:.2f}",
            color="#E74C3C",
            fontsize=12,
        )

        plt.title(f"卫星：{data['prn']} 误差变化曲线")
        plt.xlabel("时刻（秒）")
        plt.ylabel("误差（毫米）")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{data['prn']}.png"))
        plt.close()

    # 绘制广播星历拟合中误差柱状图
    plt.figure(figsize=(9, 4))
    plt.bar(prn_list, sigma_0_list, label="$\sigma_0$", color="#5DADE2")
    plt.plot(
        (-1, len(prn_list)),
        (np.mean(sigma_0_list), np.mean(sigma_0_list)),
        c="#E74C3C",
    )
    plt.text(
        -0.5,
        np.mean(sigma_0_list) + 0.5,
        f"{np.mean(sigma_0_list):.2f}",
        color="#E74C3C",
        fontsize=12,
    )

    plt.title("广播星历拟合中误差$\sigma_0$")
    plt.xlabel("卫星编号")
    plt.ylabel("误差（毫米）")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "sigma_0.png"))
    plt.close()

    # 绘制X方向的广播星历拟合误差柱状图
    plt.figure(figsize=(9, 4))
    plt.bar(prn_list, rms_x_list, label="RMSE - X", color="#2ECC71")
    plt.plot(
        (-1, len(prn_list)), (np.mean(rms_x_list), np.mean(rms_x_list)), c="#E74C3C",
    )
    plt.text(
        -0.5, np.mean(rms_x_list) + 0.5, f"{np.mean(rms_x_list):.2f}", fontsize=12,
    )

    plt.title("X方向广播星历拟合误差")
    plt.xlabel("卫星编号")
    plt.ylabel("误差（毫米）")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "rms_x.png"))
    plt.close()

    # 绘制Y方向的广播星历拟合误差柱状图
    plt.figure(figsize=(9, 4))
    plt.bar(prn_list, rms_y_list, label="RMSE - Y", color="#2980B9")
    plt.plot(
        (-1, len(prn_list)), (np.mean(rms_y_list), np.mean(rms_y_list)), c="#E74C3C",
    )
    plt.text(
        -0.5, np.mean(rms_y_list) + 0.5, f"{np.mean(rms_y_list):.2f}", fontsize=12,
    )

    plt.title("Y方向广播星历拟合误差")
    plt.xlabel("卫星编号")
    plt.ylabel("误差（毫米）")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "rms_y.png"))
    plt.close()

    # 绘制Z方向的广播星历拟合误差柱状图
    plt.figure(figsize=(9, 4))
    plt.bar(prn_list, rms_z_list, label="RMSE - Z", color="#8E44AD")
    plt.plot(
        (-1, len(prn_list)), (np.mean(rms_z_list), np.mean(rms_z_list)), c="#E74C3C",
    )
    plt.text(
        -0.5, np.mean(rms_z_list) + 0.5, f"{np.mean(rms_z_list):.2f}", fontsize=12,
    )

    plt.title("Z方向广播星历拟合误差")
    plt.xlabel("卫星编号")
    plt.ylabel("误差（毫米）")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "rms_z.png"))
    plt.close()


if __name__ == "__main__":
    main()
